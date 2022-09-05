import tornado
from server import WebSocketServer, my_function_2
from client import WebSocketClient
from datetime import datetime


def main():
    # Create a web app whose only endpoint is a WebSocket, and start the web
    # app on port 8888.
    tornado_app = tornado.web.Application(
        [(r"/websocket/", WebSocketServer)],
        websocket_ping_interval=10,
        websocket_ping_timeout=30,
    )
    tornado_app.listen(8888)

    # Create an event loop (what Tornado calls an IOLoop).
    io_loop = tornado.ioloop.IOLoop.current()

    periodic_callback = tornado.ioloop.PeriodicCallback(
        lambda: my_function_2(datetime.now().second),
        1
    )
    periodic_callback.start()

    client = WebSocketClient(io_loop)
    io_loop.add_callback(client.start)

    # Start the event loop.
    io_loop.start()


if __name__ == "__main__":
    main()
