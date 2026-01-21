#!/usr/bin/env python

# Simplified MBTiles utility for non-geographic imagery
# Based on MBUtil (c) Development Seed 2012

import sqlite3, logging, time, os, re
from datetime import datetime, timezone

logger = logging.getLogger(__name__)

def mbtiles_setup(cur):
    cur.execute("""
        create table tiles (
            plane_id integer,
            zoom_level integer,
            tile_column integer,
            tile_row integer,
            tile_data blob);
            """)
    cur.execute("""create table metadata
        (name text, value text);""")
    cur.execute("""create unique index name on metadata (name);""")
    cur.execute("""create unique index tile_index on tiles
        (plane_id, zoom_level, tile_column, tile_row);""")

def mbtiles_connect(mbtiles_file, silent):
    try:
        con = sqlite3.connect(mbtiles_file)
        return con
    except Exception as e:
        if not silent:
            logger.error("Could not connect to database")
            logger.exception(e)
        raise

def optimize_connection(cur):
    cur.execute("""PRAGMA synchronous=1""") # Set to Normal. It is safer than Off. It may slow down insertions but will prevent corrupt db
    cur.execute("""PRAGMA locking_mode=EXCLUSIVE""")
    cur.execute("""PRAGMA journal_mode=DELETE""")

def compression_prepare(cur, silent):
    if not silent:
        logger.debug('Prepare database compression.')
    cur.execute("""
      CREATE TABLE if not exists images (
        tile_data blob,
        tile_id integer);
    """)
    cur.execute("""
      CREATE TABLE if not exists map (
        plane_id integer,  
        zoom_level integer,
        tile_column integer,
        tile_row integer,
        tile_id integer);
    """)

def optimize_database(con, silent):
    if not silent:
        logger.debug('analyzing db')
    con.execute("""ANALYZE;""")
    if not silent:
        logger.debug('cleaning db')

    # VACUUM requires autocommit mode
    con.isolation_level = None
    con.execute("""VACUUM;""")
    con.isolation_level = ''


def compression_do(cur, con, chunk, silent):
    if not silent:
        logger.debug('Making database compression.')
    overlapping = 0
    unique = 0
    total = 0
    cur.execute("select count(zoom_level) from tiles")
    res = cur.fetchone()
    total_tiles = res[0]
    last_id = 0
    if not silent:
        logging.debug("%d total tiles to fetch" % total_tiles)
    for i in range(total_tiles // chunk + 1):
        if not silent:
            logging.debug("%d / %d rounds done" % (i, (total_tiles / chunk)))
        ids = []
        files = []
        start = time.time()
        cur.execute("""select plane_id, zoom_level, tile_column, tile_row, tile_data
            from tiles where rowid > ? and rowid <= ?""", ((i * chunk), ((i + 1) * chunk)))
        if not silent:
            logger.debug("select: %s" % (time.time() - start))
        rows = cur.fetchall()
        for r in rows:
            total = total + 1
            if r[4] in files:
                overlapping = overlapping + 1
                start = time.time()
                query = """insert into map
                    (plane_id, zoom_level, tile_column, tile_row, tile_id)
                    values (?, ?, ?, ?, ?)"""
                if not silent:
                    logger.debug("insert: %s" % (time.time() - start))
                # Use tile_data (r[4]) to find the existing image id
                cur.execute(query, (r[0], r[1], r[2], r[3], ids[files.index(r[4])]))
            else:
                unique = unique + 1
                last_id += 1

                ids.append(last_id)
                files.append(r[4])

                start = time.time()
                query = """insert into images
                    (tile_id, tile_data)
                    values (?, ?)"""
                cur.execute(query, (str(last_id), sqlite3.Binary(r[4])))
                if not silent:
                    logger.debug("insert into images: %s" % (time.time() - start))
                start = time.time()
                query = """insert into map
                    (plane_id, zoom_level, tile_column, tile_row, tile_id)
                    values (?, ?, ?, ?, ?)"""
                cur.execute(query, (r[0], r[1], r[2], r[3], last_id))
                if not silent:
                    logger.debug("insert into map: %s" % (time.time() - start))
        con.commit()

    if not silent:
        dedup_pct = (overlapping / total * 100) if total > 0 else 0
        logger.info('Compression: %d tiles processed, %d unique, %d duplicates (%.1f%% deduplication)' %
                   (total, unique, overlapping, dedup_pct))

def compression_finalize(cur, con, silent):
    if not silent:
        logger.debug('Finalizing database compression.')
    cur.execute("""drop table tiles;""")
    cur.execute("""create view tiles as
        select map.plane_id as plane_id, 
               map.zoom_level as zoom_level,
               map.tile_column as tile_column,
               map.tile_row as tile_row,
               images.tile_data as tile_data FROM 
            map JOIN images on images.tile_id = map.tile_id;""")
    cur.execute("""
          CREATE UNIQUE INDEX map_index on map
            (plane_id, zoom_level, tile_column, tile_row);""")
    cur.execute("""
          CREATE UNIQUE INDEX images_id on images
            (tile_id);""")

    # VACUUM requires autocommit mode
    con.isolation_level = None
    cur.execute("""vacuum;""")
    con.isolation_level = ''

    cur.execute("""analyze;""")

def get_dirs(path):
    return [name for name in os.listdir(path)
        if os.path.isdir(os.path.join(path, name))]

def _parse_plane_id(name: str, silent: bool) -> int:
    """Extract an integer plane_id from a directory name.
    Accepts names like '123', 'd_123', 'plane123'. Uses trailing digits.
    Raises ValueError if no digits are found.
    """
    name = name.strip()
    # if the whole name is digits
    if name.isdigit():
        return int(name)
    # find trailing digits (e.g., d_123, plane123)
    m = re.search(r"(\d+)$", name)
    if m:
        return int(m.group(1))
    if not silent:
        logger.error("Cannot parse plane_id from directory '%s'. Expected trailing integer (e.g., 'd_123' or '123').", name)
    raise ValueError(f"Invalid plane directory name: {name}")

def disk_to_mbtiles(directory_path, mbtiles_file, **kwargs):
    """
    Import tiles from disk into MBTiles database.
    Expects tiles in plane/z/y/x.jpg layout (non-geographic imagery).

    Args:
        directory_path: Path to tile directory
        mbtiles_file: Output MBTiles file path
        **kwargs:
            format: Tile format (default: 'png')
            batch_size: Tiles per batch insert (default: 1000)
            compression: Enable tile deduplication (default: False)
            compression_chunk: Tiles per compression round (default: 10000)
            silent: Suppress logging (default: False)
            name: Dataset name (optional)
            description: Dataset description (optional)
            width: Image width in pixels (optional)
            height: Image height in pixels (optional)
    """
    silent = kwargs.get('silent', False)
    image_format = kwargs.get('format', 'png')
    image_format = (image_format or 'png').lower()
    batch_size = kwargs.get('batch_size', 1000)
    sample_every = max(10000, batch_size * 10)

    if not silent:
        logger.info("Importing disk to MBTiles")
        logger.debug("%s --> %s" % (directory_path, mbtiles_file))

    con = mbtiles_connect(mbtiles_file, silent)
    cur = con.cursor()
    optimize_connection(cur)
    mbtiles_setup(cur)

    count = 0
    start_time = time.time()
    batch = []

    # Track metadata while iterating
    plane_ids = set()
    zoom_levels = set()

    # Iterate through plane/z/y/x directory structure
    for plane_dir in get_dirs(directory_path):
        try:
            d = _parse_plane_id(plane_dir, silent)
        except ValueError:
            if not silent:
                logger.warning("Skipping directory without numeric plane id: %s", plane_dir)
            continue
        plane_ids.add(d)
        for zoom_dir in get_dirs(os.path.join(directory_path, plane_dir)):
            z = int(zoom_dir)
            zoom_levels.add(z)
            for row_dir in get_dirs(os.path.join(directory_path, plane_dir, zoom_dir)):
                y = int(row_dir)  # y coordinate from directory name
                for current_file in os.listdir(os.path.join(directory_path, plane_dir, zoom_dir, row_dir)):
                    if current_file == ".DS_Store":
                        if not silent:
                            logger.debug("Skipping .DS_Store file")
                        continue

                    # Robust extension parsing (handles multi-dot names, case-insensitive)
                    file_name, ext = os.path.splitext(current_file)
                    ext = ext.lstrip('.').lower()
                    if ext != image_format or not file_name:
                        continue

                    x = int(file_name)  # x coordinate from file name

                    file_path = os.path.join(directory_path, plane_dir, zoom_dir, row_dir, current_file)
                    with open(file_path, 'rb') as f:
                        file_content = f.read()

                    # Log only the first few and then sparsely to avoid huge logs
                    if not silent and (count < 10 or ((count + 1) % sample_every == 0)):
                        logger.debug(' Read tile from Plane (d): %i, Zoom (z): %i\tCol (x): %i\tRow (y): %i' % (d, z, x, y))

                    # Add to batch
                    batch.append((d, z, x, y, sqlite3.Binary(file_content)))
                    count += 1

                    # Insert batch when it reaches batch_size
                    if len(batch) >= batch_size:
                        cur.executemany("""insert into tiles (
                            plane_id, zoom_level, tile_column, tile_row, tile_data) values
                            (?, ?, ?, ?, ?);""", batch)
                        con.commit()
                        batch = []
                        if not silent:
                            logger.info(" %s tiles inserted (%d tiles/sec)" % (count, count / (time.time() - start_time)))

    # Insert remaining tiles in batch
    if batch:
        cur.executemany("""insert into tiles (
            plane_id, zoom_level, tile_column, tile_row, tile_data) values
            (?, ?, ?, ?, ?);""", batch)
        con.commit()

    if not silent:
        logger.info('Total: %s tiles inserted in %.2f seconds (%d tiles/sec)' %
                   (count, time.time() - start_time, count / (time.time() - start_time)))

    # Insert metadata
    metadata = {
        # Auto-detected
        'format': image_format,
        'minzoom': str(min(zoom_levels)) if zoom_levels else '0',
        'maxzoom': str(max(zoom_levels)) if zoom_levels else '0',
        'planes': ','.join(str(p) for p in sorted(plane_ids)),
        'plane_count': str(len(plane_ids)),
        'created': datetime.now(timezone.utc).isoformat(),
        'tile_count': str(count),
    }
    # Optional (from kwargs)
    if kwargs.get('name'):
        metadata['name'] = kwargs['name']
    if kwargs.get('description'):
        metadata['description'] = kwargs['description']
    if kwargs.get('width'):
        metadata['width'] = str(kwargs['width'])
    if kwargs.get('height'):
        metadata['height'] = str(kwargs['height'])
    if kwargs.get('voxel_size'):
        # Store voxel size as comma-separated string: "x,y,z"
        voxel_size = kwargs['voxel_size']
        if isinstance(voxel_size, (list, tuple)) and len(voxel_size) == 3:
            metadata['voxel_size'] = ','.join(str(v) for v in voxel_size)

    for name, value in metadata.items():
        cur.execute('INSERT INTO metadata (name, value) VALUES (?, ?)', (name, value))
    con.commit()

    if not silent:
        logger.info('Metadata: %s', metadata)

    if kwargs.get('compression', False):
        compression_prepare(cur, silent)
        compression_chunk = kwargs.get('compression_chunk', 10000)
        compression_do(cur, con, compression_chunk, silent)
        compression_finalize(cur, con, silent)

    optimize_database(con, silent)

    con.close()
