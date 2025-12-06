#!/usr/bin/env python

# Simplified MBTiles utility for non-geographic imagery
# Based on MBUtil (c) Development Seed 2012

import sqlite3, logging, time, os, re
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
from queue import Queue
from threading import Thread

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
    cur.execute("""PRAGMA synchronous=0""")
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


# =============================================================================
# PARALLEL VERSION - Producer/Consumer Pattern
# =============================================================================
#
# Architecture:
#   PRODUCERS (N threads)          QUEUE              CONSUMER (1 thread)
#   ┌─────────────────┐           ┌─────┐           ┌─────────────────┐
#   │ Read tile files │ ───────►  │     │  ───────► │ Batch insert    │
#   │ from disk       │           │     │           │ into SQLite     │
#   └─────────────────┘           └─────┘           └─────────────────┘
#
# Why this works:
#   - File reading benefits from parallelism (multiple I/O operations)
#   - SQLite writing stays single-threaded (avoids lock contention)
#   - Queue provides backpressure if writers fall behind readers
# =============================================================================


def _collect_tile_paths(directory_path, image_format, silent):
    """
    Walk directory tree and collect all tile file paths.
    This is fast (no file I/O, just directory listing).

    Returns:
        List of (plane_id, z, x, y, file_path) tuples
    """
    tile_paths = []
    for plane_dir in get_dirs(directory_path):
        try:
            plane_id = _parse_plane_id(plane_dir, silent)
        except ValueError:
            if not silent:
                logger.warning("Skipping directory: %s", plane_dir)
            continue

        plane_path = os.path.join(directory_path, plane_dir)
        for zoom_dir in get_dirs(plane_path):
            z = int(zoom_dir)
            zoom_path = os.path.join(plane_path, zoom_dir)

            for row_dir in get_dirs(zoom_path):
                y = int(row_dir)
                row_path = os.path.join(zoom_path, row_dir)

                for filename in os.listdir(row_path):
                    if filename == ".DS_Store":
                        continue
                    name, ext = os.path.splitext(filename)
                    if ext.lstrip('.').lower() != image_format or not name:
                        continue

                    x = int(name)
                    file_path = os.path.join(row_path, filename)
                    tile_paths.append((plane_id, z, x, y, file_path))

    return tile_paths


def _read_tile(tile_info):
    """
    PRODUCER: Read a single tile file from disk.

    Args:
        tile_info: (plane_id, z, x, y, file_path)

    Returns:
        (plane_id, z, x, y, file_bytes)
    """
    plane_id, z, x, y, file_path = tile_info
    with open(file_path, 'rb') as f:
        data = f.read()
    return (plane_id, z, x, y, data)


def _writer_worker(queue, mbtiles_file, batch_size, silent):
    """
    CONSUMER: Take tiles from queue and batch-insert into SQLite.
    Runs in a dedicated thread. Stops when it receives None (poison pill).

    Args:
        queue: Queue to read from
        mbtiles_file: Path to MBTiles database
        batch_size: Number of tiles per batch insert
        silent: Suppress logging
    """
    con = sqlite3.connect(mbtiles_file)
    cur = con.cursor()
    optimize_connection(cur)

    batch = []
    count = 0
    start_time = time.time()

    while True:
        item = queue.get()

        # Poison pill signals shutdown
        if item is None:
            break

        plane_id, z, x, y, data = item
        batch.append((plane_id, z, x, y, sqlite3.Binary(data)))
        count += 1

        if len(batch) >= batch_size:
            cur.executemany(
                """INSERT INTO tiles (plane_id, zoom_level, tile_column, tile_row, tile_data)
                   VALUES (?, ?, ?, ?, ?)""",
                batch
            )
            con.commit()
            if not silent:
                elapsed = time.time() - start_time
                logger.info(" %d tiles written (%.0f tiles/sec)" % (count, count / elapsed))
            batch = []

    # Insert remaining tiles
    if batch:
        cur.executemany(
            """INSERT INTO tiles (plane_id, zoom_level, tile_column, tile_row, tile_data)
               VALUES (?, ?, ?, ?, ?)""",
            batch
        )
        con.commit()

    con.close()
    return count


def disk_to_mbtiles_parallel(directory_path, mbtiles_file, **kwargs):
    """
    Import tiles from disk into MBTiles using parallel file reading.

    Same interface as disk_to_mbtiles, but uses producer-consumer pattern:
    - Multiple threads read tile files from disk (producers)
    - Single thread writes to SQLite database (consumer)

    Args:
        directory_path: Path to tile directory
        mbtiles_file: Output MBTiles file path
        **kwargs:
            format: Tile format (default: 'png')
            batch_size: Tiles per batch insert (default: 50000)
            workers: Number of reader threads (default: 8)
            silent: Suppress logging (default: False)
            name: Dataset name (optional)
            description: Dataset description (optional)
            width: Image width in pixels (optional)
            height: Image height in pixels (optional)
    """
    silent = kwargs.get('silent', False)
    image_format = kwargs.get('format', 'png').lower()
    batch_size = kwargs.get('batch_size', 50000)
    workers = kwargs.get('workers', 8)

    if not silent:
        logger.info("Importing disk to MBTiles (parallel: %d workers)" % workers)

    # -------------------------------------------------------------------------
    # Step 1: Collect all tile paths (fast, no file I/O)
    # -------------------------------------------------------------------------
    if not silent:
        logger.info("Step 1/4: Collecting tile paths...")
    tile_paths = _collect_tile_paths(directory_path, image_format, silent)
    total_tiles = len(tile_paths)

    if not silent:
        logger.info("Found %d tiles to process" % total_tiles)

    if total_tiles == 0:
        logger.warning("No tiles found!")
        return

    # Extract metadata from paths
    plane_ids = set(t[0] for t in tile_paths)
    zoom_levels = set(t[1] for t in tile_paths)

    # -------------------------------------------------------------------------
    # Step 2: Setup database
    # -------------------------------------------------------------------------
    if not silent:
        logger.info("Step 2/4: Setting up database...")
    con = mbtiles_connect(mbtiles_file, silent)
    cur = con.cursor()
    optimize_connection(cur)
    mbtiles_setup(cur)
    con.commit()
    con.close()  # Close - writer thread will open its own connection

    # -------------------------------------------------------------------------
    # Step 3: Start consumer thread, then produce tiles
    # -------------------------------------------------------------------------
    if not silent:
        logger.info("Step 3/4: Reading and writing tiles...")

    queue = Queue(maxsize=batch_size * 2)  # Backpressure if consumer falls behind

    # Start consumer thread
    writer = Thread(
        target=_writer_worker,
        args=(queue, mbtiles_file, batch_size, silent),
        name="sqlite-writer"
    )
    writer.start()

    # Produce tiles using thread pool
    start_time = time.time()
    with ThreadPoolExecutor(max_workers=workers) as executor:
        for result in executor.map(_read_tile, tile_paths):
            queue.put(result)

    # Signal consumer to stop
    queue.put(None)
    writer.join()

    elapsed = time.time() - start_time
    if not silent:
        logger.info("Total: %d tiles in %.2f seconds (%.0f tiles/sec)" %
                   (total_tiles, elapsed, total_tiles / elapsed))

    # -------------------------------------------------------------------------
    # Step 4: Add metadata
    # -------------------------------------------------------------------------
    if not silent:
        logger.info("Step 4/4: Writing metadata...")

    con = mbtiles_connect(mbtiles_file, silent)
    cur = con.cursor()

    metadata = {
        'format': image_format,
        'minzoom': str(min(zoom_levels)) if zoom_levels else '0',
        'maxzoom': str(max(zoom_levels)) if zoom_levels else '0',
        'planes': ','.join(str(p) for p in sorted(plane_ids)),
        'plane_count': str(len(plane_ids)),
        'created': datetime.now(timezone.utc).isoformat(),
        'tile_count': str(total_tiles),
    }
    if kwargs.get('name'):
        metadata['name'] = kwargs['name']
    if kwargs.get('description'):
        metadata['description'] = kwargs['description']
    if kwargs.get('width'):
        metadata['width'] = str(kwargs['width'])
    if kwargs.get('height'):
        metadata['height'] = str(kwargs['height'])

    for name, value in metadata.items():
        cur.execute('INSERT INTO metadata (name, value) VALUES (?, ?)', (name, value))
    con.commit()

    if not silent:
        logger.info('Metadata: %s', metadata)

    optimize_database(con, silent)
    con.close()
