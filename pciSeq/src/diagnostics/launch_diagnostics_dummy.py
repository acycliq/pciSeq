import os
import asyncio
import streamlit.bootstrap
from streamlit.server.server import Server
import tornado.ioloop
import sqlite3
import pandas as pd
from multiprocessing import Process, Queue
from datetime import datetime
from tornado.platform.asyncio import AsyncIOLoop
from streamlit import config as _config
import nest_asyncio

# nest_asyncio.apply()

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, 'diagnostics_dummy.py')

# _config.set_option("server.headless", True)
args = []



async def do_some_iterations(con):
    # con = sqlite3.connect("file:memdb1?mode=memory&cache=shared", uri=True)
    # df = pd.DataFrame({'x': [11, 21, 31], 'y': ['a', 'b', 'c']})
    # df.to_sql(name='spots', con=con, if_exists='replace', index=False)
    # con.close()
    # print('Database created')
    for i in range(10):
        print(datetime.now().time())
        df = pd.DataFrame({'run': [i]})
        df.to_sql(name='spots', con=con, if_exists='replace', index=False)
        await asyncio.sleep(1)
    print('... Cool!')


async def main():
    con = sqlite3.connect("file:memdb1?mode=memory&cache=shared")
    await asyncio.gather(
        do_some_iterations(con),
        do_something(con)
    )


def run(
    main_script_path: str,
        command_line,
        args,
        flag_options
) -> None:
    streamlit.bootstrap._fix_sys_path(main_script_path)
    streamlit.bootstrap._fix_matplotlib_crash()
    streamlit.bootstrap._fix_tornado_crash()
    streamlit.bootstrap._fix_sys_argv(main_script_path, args)
    streamlit.bootstrap._fix_pydeck_mapbox_api_warning()
    streamlit.bootstrap._install_config_watchers(flag_options)
    streamlit.bootstrap._install_pages_watcher(main_script_path)

    # Install a signal handler that will shut down the ioloop
    # and close all our threads
    streamlit.bootstrap._set_up_signal_handler()

    ioloop = tornado.ioloop.IOLoop.current()

    # Create and start the server.
    server = Server(ioloop, main_script_path, command_line)
    server.start(streamlit.bootstrap._on_server_start)

    print('im here')
    # Start the ioloop. This function will not return until the
    # server is shut down.
    # ioloop.start()


async def do_something(con):
    print('do...')
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, 'diagnostics_dummy.py')

    # _config.set_option("server.headless", True)
    args = []
    run(filename, '', args, flag_options={})
    print('...somethinggg')


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except:
        pass
    # do_something()

    # loop = asyncio.get_event_loop()
    # # Blocking call which returns when the display_date() coroutine is done
    # loop.run_until_complete(do_some_iterations())
    # print('End')
    # loop.close()
    # print('End')

# streamlit.bootstrap.run(filename, '', args, flag_options={})