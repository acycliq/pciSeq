import os
import streamlit.bootstrap
from streamlit.server.server import Server
import tornado.ioloop


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


async def launch_dashboard():
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, 'diagnostics.py')

    # _config.set_option("server.headless", True)
    args = []
    run(filename, '', args, flag_options={})


if __name__ == "__main__":
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, 'diagnostics.py')

    # _config.set_option("server.headless", True)
    args = []
    streamlit.bootstrap.run(filename, '', args, flag_options={})
