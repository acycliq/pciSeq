""" Stream data from the WebSocket and update the Beta posterior parameters online. """

import tornado.ioloop
import tornado.websocket


class WebSocketClient:
    def __init__(self, io_loop):
        self.connection = None
        self.io_loop = io_loop
        self.num_successes = 0
        self.num_trials = 0

    def start(self):
        self.connect_and_read()

    def stop(self):
        self.io_loop.stop()

    def connect_and_read(self):
        print("Reading...")
        tornado.websocket.websocket_connect(
            url=f"ws://localhost:8888/websocket/",
            callback=self.maybe_retry_connection,
            on_message_callback=self.on_message,
            ping_interval=10,
            ping_timeout=30,
        )

    def maybe_retry_connection(self, future) -> None:
        try:
            self.connection = future.result()
        except:
            print("Could not reconnect, retrying in 3 seconds...")
            self.io_loop.call_later(3, self.connect_and_read)

    def on_message(self, message):
        if message is None:
            print("Disconnected, reconnecting...")
            self.connect_and_read()

        # print(message)
        message = int(message)
        # print('message is %3.2f' % message)
        # message = int('12')
        self.num_successes += message
        self.num_trials += 1

        alpha = 2 + self.num_successes
        beta = 2 + self.num_trials - self.num_successes
        mean = self.num_successes / self.num_trials
        print("alpha = %d, beta = %d, mean = %3.2f" % (alpha, beta, mean))
        # print(f"α = {alpha}; β = {beta}; mean = {mean}")


def main():
    # Create an event loop (what Tornado calls an IOLoop).
    io_loop = tornado.ioloop.IOLoop.current()

    # Before starting the event loop, instantiate a WebSocketClient and add a
    # callback to the event loop to start it. This way the first thing the
    # event loop does is to start the client.
    client = WebSocketClient(io_loop)
    io_loop.add_callback(client.start)

    # Start the event loop.
    io_loop.start()


if __name__ == "__main__":
    main()