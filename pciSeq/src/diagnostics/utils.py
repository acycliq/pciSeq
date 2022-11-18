import os
import redis
import pickle
from sys import platform
from pciSeq.src.cell_call.log_config import attach_to_log, logger
import pciSeq.src.diagnostics.config as config
import getpass
import subprocess as sp
import sysconfig


class redis_db():
    def __init__(self, flush=False):
        self.pool = None
        self.redis_client = None
        self.is_connected = None
        self._attach_client()
        if flush:
            self._flushdb()

    def _attach_client(self):
        self.pool = config.SETTINGS['POOL']
        self.redis_client = redis.Redis(connection_pool=self.pool)
        self.is_connected = self._is_connected()

    def to_redis(self, df_in, key, **kwargs):
        df = df_in.copy()
        df = df.assign(**kwargs)
        self.redis_client.set(key, pickle.dumps(df))

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
            code_1, code_2 = validate()
            if code_1 == 0 and code_2 == 0:
                return True
            else:
                raise Exception("Cannot validate redis")

    # def is_connected(self):
    #     try:
    #         return self.redis_client.ping()
    #     except (redis.exceptions.ConnectionError, ConnectionRefusedError):
    #         if check_platform() == "windows" and confirm_prompt("Do you want to instal Memurai?"):
    #             logger.info("Redis ping failed!. Trying to install redis server")
    #             os.path.join(os.path.dirname(os.path.realpath(__file__)), )
    #             msi = os.path.join(sysconfig.get_path('purelib'), 'pciSeq', 'static', 'memurai',
    #                                'Memurai-Developer-v3.1.4.msi')
    #             if not os.path.isfile(msi):
    #                 msi = os.path.join(os.getcwd(), 'static', 'memurai', 'Memurai-Developer-v3.1.4.msi')
    #             logger.info("Calling %s" % msi)
    #             exit_code = sp.call(msi, shell=True)  # returns the exit code in unix. if 0 then success
    #             if exit_code > 0:
    #                 raise Exception("Installation failed with exit code: %d" % exit_code)
    #             return exit_code == 0
    #         elif check_platform() in ["linux", "osx"] and confirm_prompt(
    #                 "Do you want to instal redis-server and redis-tools?"):
    #             sudo_password = getpass.getpass(prompt='sudo password: ')
    #             subprocess_cmd(['sudo', '-S', 'apt-get', 'install', '-y', 'redis-server', 'redis-tools'], sudo_password)
    #             subprocess_cmd(['sudo', 'service', 'redis-server', 'start'], sudo_password)
    #             logger.info("Pinging redis: ...")
    #             out, err = subprocess_cmd(['redis-cli', 'ping'], sudo_password)
    #             if out.decode('UTF-8').rstrip() == 'PONG':
    #                 logger.info("%s" % out.decode('UTF-8'))
    #             else:
    #                 raise Exception("%s" % out)
    #             return True
    #     except Exception as e:
    #         raise Exception


def is_running(os):
    out = None
    err = None
    exit_code = None
    if os == "windows":
        out, err, exit_code = subprocess_cmd(['memurai-cli.exe', 'ping'])
    elif os in ['linux', 'osx']:
        out, err, exit_code = subprocess_cmd(['redis-cli', 'ping'])
    else:
        raise Exception('not implemented')
    if not err.decode('UTF-8') == '':
        logger.info(err.decode('UTF-8'))
        if confirm_prompt("Server is not running, do you want to start it?"):
            out, err, exit_code = start_server(os)
    return out, err, exit_code


def start_server(os):
    out = None
    err = None
    exit_code = None
    if os == "windows":
        out, err, exit_code = subprocess_cmd(['memurai.exe', '--service-start'])
    elif os in ['linux', 'osx']:
        passwd = getpass.getpass(prompt='sudo password: ')
        out, err, exit_code = subprocess_cmd(['sudo', '-S', 'service', 'redis-server', 'start'], passwd)
    else:
        raise Exception('not implemented')
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
    if os == "windows":
        out, err, exit_code = subprocess_cmd(['memurai.exe', '--version'])
    elif os in ['linux', 'osx']:
        out, err, exit_code = subprocess_cmd(['redis-server', '--version'])
    else:
        raise Exception('not implemented')
    if not err.decode('UTF-8') == '':
        logger.info(err.decode('UTF-8').rstrip())
        if confirm_prompt("Server is not installed, do you want to install it?"):
            install_redis_server(os)
    return out, err, exit_code


def install_redis_server(os):
    out = None
    err = None
    exit_code = None
    if os == "windows":
        msi = locate_msi()
        logger.info("Calling %s" % msi)
        out, err, exit_code = subprocess_cmd([msi])
    elif os in ['linux', 'osx']:
        sudo_password = getpass.getpass(prompt='sudo password: ')
        out, err, exit_code = subprocess_cmd(['sudo', '-S', 'apt-get', 'install', '-y', 'redis-server', 'redis-tools'], sudo_password)
    else:
        raise Exception('not implemented')
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
    msi = os.path.join(sysconfig.get_path('purelib'), 'pciSeq', 'static', 'memurai', 'Memurai-Developer-v3.1.4.msi')
    if not os.path.isfile(msi):
        msi = os.path.join(os.getcwd(), 'static', 'memurai', 'Memurai-Developer-v3.1.4.msi')
    return msi


def confirm_prompt(question):
    reply = None
    while reply not in ("", "y", "n"):
        reply = input(f"{question} (y/n): ").lower()
    return (reply in ("", "y"))


def check_platform():
    if platform == "linux" or platform == "linux2":
        return "linux"
    elif platform == "darwin":
        return "osx"
    elif platform == "win32":
        return "windows"
    else:
        raise Exception("Cannot determine OS")


def validate():
    os = check_platform()
    _, _,  code_1 = is_redis_installed(os)
    _, _,  code_2 = is_running(os)
    return code_1, code_2


if __name__ == "__main__":
    attach_to_log()
    logger.info("hwllo")
    # os = "windows"
    os = "linux"
    is_redis_installed(os)
    is_running(os)
    stop_server(os)
    print('done!')

