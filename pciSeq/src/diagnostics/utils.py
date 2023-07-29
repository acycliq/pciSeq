import os
import redis
import pickle
from sys import platform
from pciSeq.src.core.utils import get_pciSeq_install_dir
from pciSeq.src.core.log_config import attach_to_log, logger
import pciSeq.src.diagnostics.config as config
import getpass
import subprocess as sp


class redis_db():
    def __init__(self, flush=False):
        self.pool = None
        self.redis_client = None
        self.is_connected = None
        self.keyspace_events_enabled = None
        self._attach_client()
        if flush:
            self._flushdb()

    def _attach_client(self):
        self.pool = config.SETTINGS['POOL']
        self.redis_client = redis.Redis(connection_pool=self.pool)
        self.is_connected = self._is_connected()
        if self.is_connected and not self.keyspace_events_enabled:
            self.enable_keyspace_events()

    def to_redis(self, df_in, key, **kwargs):
        df = df_in.copy()
        df = df.assign(**kwargs)
        self.redis_client.set(key, pickle.dumps(df))

    def publish(self, df_in, channel, **kwargs):
        df = df_in.copy()
        df = df.assign(**kwargs)
        self.redis_client.publish(channel, pickle.dumps(df))

    def from_redis(self, key):
        """Retrieve Numpy array from Redis key 'key'"""
        return pickle.loads(self.redis_client.get(key))

    def get_db_tables(self):
        bytestings = self.redis_client.keys()
        return [d.decode('UTF-8') for d in bytestings]

    def _flushdb(self):
        self.redis_client.flushdb()

    def _is_connected(self):
        try:
            return self.redis_client.ping()
        except (redis.exceptions.ConnectionError, ConnectionRefusedError):
            code_1, code_2 = validate_redis()
            if code_1 == 0 and code_2 == 0:
                return True
            else:
                raise Exception("Cannot validate redis")

    def enable_keyspace_events(self):
        exe = "memurai-cli.exe" if check_platform() == "windows" else "redis-cli"
        try:
            out, err, exit_code = subprocess_cmd([exe, 'config', 'set', 'notify-keyspace-events', 'KEA'])
            if exit_code != 0:
                logger.info(out.decode('UTF-8').rstrip())
                logger.info(err.decode('UTF-8').rstrip())
                raise Exception('notify-keyspace-events failed with exit code: %d' % exit_code)
            logger.info(" enabling keyspace events... %s" % out.decode('UTF-8').rstrip())
            self.keyspace_events_enabled = True
        except OSError as ex:
            logger.info("Cannot enable keyspace events. Failed with error: %s" % ex)
            raise


def is_redis_running(os):
    out = None
    err = None
    exit_code = None
    exe = "memurai-cli.exe" if check_platform() == "windows" else "redis-cli"
    out, err, exit_code = subprocess_cmd([exe, 'ping'])
    if exit_code != 0:
        logger.info(err.decode('UTF-8'))
        logger.info(" Starting redis server ...")
        out, err, exit_code = start_server(os)
    return out, err, exit_code


def start_server(os):
    out = None
    err = None
    exit_code = None
    if os == "windows":
        out, err, exit_code = subprocess_cmd(['memurai.exe', '--service-start'])
    elif os in ['linux']:
        passwd = getpass.getpass(prompt='sudo password: ')
        out, err, exit_code = subprocess_cmd(['sudo', '-S', 'service', 'redis-server', 'start'], passwd)
    else:
        ## need to add osx here!
        raise Exception('Cannot automatically start redis server under %s. Not implemented. Try manually.' % os)
    if not out.decode('UTF-8') == '':
        logger.info(out.decode('UTF-8').rstrip())
    return out, err, exit_code


def stop_server(os):
    out = None
    err = None
    exit_code = None
    if os == "windows":
        out, err, exit_code = subprocess_cmd(['memurai.exe', '--service-stop'])
    elif os in ['linux', 'osx']:
        passwd = getpass.getpass(prompt='sudo password: ')
        out, err, exit_code = subprocess_cmd(['sudo', '-S', 'service', 'redis-server', 'stop'], passwd)
    else:
        raise Exception('not implemented')
    if not out.decode('UTF-8') == '':
        logger.info(out.decode('UTF-8').rstrip())
    return out, err, exit_code


def is_redis_installed(os):
    out = None
    err = None
    exit_code = None

    exe = "memurai.exe" if check_platform() == "windows" else "redis-server"
    out, err, exit_code = subprocess_cmd([exe, '--version'])
    if not err.decode('UTF-8') == '':
        logger.info(err.decode('UTF-8').rstrip())
        if confirm_prompt("Server is not installed, do you want to install it?"):
            out, err, exit_code = install_redis_server(os)
    return out, err, exit_code


def install_redis_server(os):
    out = None
    err = None
    exit_code = None
    if os == "windows":
        msi = locate_msi()
        logger.info("Calling %s" % msi)
        out, err, exit_code = subprocess_cmd([msi])
    elif os in ['linux']:
        sudo_password = getpass.getpass(prompt='sudo password: ')
        out, err, exit_code = subprocess_cmd(['sudo', '-S', 'apt-get', 'install', '-y', 'redis-server', 'redis-tools'], sudo_password)
    else:
        ## need to add osx here!
        raise Exception('not implemented for %s' % os)
    if exit_code > 0:
        raise Exception("installation failed with exit code %d" % exit_code)
    return out, err, exit_code


def subprocess_cmd(command, passwd=None):
    p = sp.Popen(command, stderr=sp.PIPE, stdout=sp.PIPE, stdin=sp.PIPE)
    out = None
    err = None
    try:
        if passwd:
            out, err = p.communicate(input=(passwd + '\n').encode(), timeout=5)
        else:
            out, err = p.communicate(timeout=5)
    except sp.TimeoutExpired:
        p.kill()
    if not err.decode('UTF-8') == '':
        p.kill()
        # logger.info(err.decode('UTF-8'))
        # raise Exception(err.decode('UTF-8'))
    return out, err, p.returncode


def locate_msi():
    pciSeq_dir = get_pciSeq_install_dir()
    msi = os.path.join(pciSeq_dir, 'static', 'memurai', 'Memurai-Developer-v3.1.4.msi')
    # msi = os.path.join(sysconfig.get_path('purelib'), 'pciSeq', 'static', 'memurai', 'Memurai-Developer-v3.1.4.msi')
    if not os.path.isfile(msi):
        msi = os.path.join(os.getcwd(), 'static', 'memurai', 'Memurai-Developer-v3.1.4.msi')
    return msi


def confirm_prompt(question):
    reply = None
    while reply not in ("", "y", "n"):
        reply = input(f"{question} (y/n): ").lower()
    return reply in ("", "y")


def check_platform():
    if platform == "linux" or platform == "linux2":
        return "linux"
    elif platform == "darwin":
        return "osx"
    elif platform == "win32":
        return "windows"
    else:
        raise Exception("Cannot determine OS")


def validate_redis():
    os = check_platform()
    _, _,  code_1 = is_redis_installed(os)
    _, _,  code_2 = is_redis_running(os)
    return code_1, code_2


if __name__ == "__main__":
    attach_to_log()
    logger.info("hwllo")
    # os = "windows"
    os = "linux"
    is_redis_installed(os)
    is_redis_running(os)
    stop_server(os)
    print('done!')

