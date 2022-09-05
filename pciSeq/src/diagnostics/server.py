""" Every 100ms, sample from a Bernoulli and write the value to a WebSocket. """

import random
import tornado.ioloop
import tornado.web
import tornado.websocket
from datetime import datetime


class WebSocketServer(tornado.websocket.WebSocketHandler):
    """Simple WebSocket handler to serve clients."""

    # Note that `clients` is a class variable and `send_message` is a
    # classmethod.
    clients = set()

    def open(self):
        # print('added client')
        WebSocketServer.clients.add(self)
        print('added client %d' % len(WebSocketServer.clients))

    def on_close(self):
        WebSocketServer.clients.remove(self)

    @classmethod
    def send_message(cls, message: str):
        print(f"Sending message {message} to {len(cls.clients)} client(s).")
        for client in cls.clients:
            client.write_message(message)


class RandomBernoulli:
    def __init__(self):
        self.p = 0.72
        print(f"True p = {self.p}")

    def sample(self):
        return int(random.uniform(0, 1) <= self.p)


def my_function(x):
    if x % 10 == 0:
        WebSocketServer.send_message(str(x))


def my_function_2(x):
    for i in range(1000000):
        if i % 100 == 1:
            print(i)
            WebSocketServer.send_message(str(i))




def main():
    # Create a web app whose only endpoint is a WebSocket, and start the web
    # app on port 8888.
    app = tornado.web.Application(
        [(r"/websocket/", WebSocketServer)],
        websocket_ping_interval=10,
        websocket_ping_timeout=30,
    )
    app.listen(8888)

    # Create an event loop (what Tornado calls an IOLoop).
    io_loop = tornado.ioloop.IOLoop.current()

    # Before starting the event loop, instantiate a RandomBernoulli and
    # register a periodic callback to write a sampled value to the WebSocket
    # every 100ms.
    # random_bernoulli = RandomBernoulli()

    # periodic_callback = tornado.ioloop.PeriodicCallback(
    #     lambda: WebSocketServer.send_message(str(random_bernoulli.sample())), 100
    # )
    # periodic_callback.start()
    #
    periodic_callback = tornado.ioloop.PeriodicCallback(
        lambda: my_function_2(datetime.now().second),
        1000
    )
    periodic_callback.start()

    # Start the event loop.
    io_loop.start()


if __name__ == "__main__":
    main()