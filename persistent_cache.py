# Implements a decorator, that caches function call values cross runtimes

from typing import Any
import atexit
import pickle
import os.path
import os
from time import sleep


print("Using persistent caching provided by persistent_cache.py")


def retry(function):
    timeout = 0.1
    retries = 10000

    def retry_wrapper(*args):
        tries = 0
        while tries < retries:
            try:
                value = function(*args)
                return value
            except:  # noqa: E722
                sleep(timeout)
                tries += 1
                if tries % 1000 == 0:
                    print(f"Function {function.__name__} is blocked!!!")


def persistent_cache(
    path_arg_mask: tuple[int, ...] = (),
    filename: str | None = None,
):
    if type(path_arg_mask) is not tuple:
        # for single argument path masks
        path_arg_mask = tuple((path_arg_mask,))

    def persistent_cache_decorator(function):
        # actual decorator with given parameters
        loaded_cache = {}

        def persistent_cache_wrapper(*args):
            path_args = path_arguments(args)
            if path_args not in loaded_cache:
                load(path_args)

            if args in loaded_cache[path_args]:
                return loaded_cache[path_args][args]
            else:
                print("[persistent_cache]: Running function")
                rv = function(*args)
                loaded_cache[path_args][args] = rv
                return rv

        def path_arguments(args: tuple[Any, ...]) -> tuple[Any, ...]:
            return tuple((args[i] for i in path_arg_mask))

        def get_path(path_args):
            p = ".persistent_cache"
            if filename is not None:
                p = os.path.join(p, filename)
            p = os.path.join(p, function.__name__)
            p = os.path.join(p, *(str(path_segment) for path_segment in path_args))
            p += ".cache"
            return p

        @retry
        def lock(lock_path):
            with open(lock_path, "x"):
                pass

        def save():
            # save files

            for path_args, sub_cache in loaded_cache.items():
                # create dir if it does not exists
                file_path = get_path(path_args)
                dir_path = os.path.dirname(file_path)
                if not os.path.isdir(dir_path):
                    os.makedirs(dir_path)
                # lock file
                lock_file_path = file_path + ".lock"
                lock(lock_file_path)
                # get newest content of cache file
                if os.path.exists(file_path):
                    with open(file_path, "rb") as cache_file:
                        combined_cache = pickle.load(cache_file)
                        # combine with new caches. If conflicts happen, the newer value is assumed to be correct
                        combined_cache.update(sub_cache)
                else:
                    combined_cache = sub_cache
                # save
                with open(file_path, "wb") as cache_file:
                    pickle.dump(combined_cache, cache_file)
                # unlock file
                os.remove(lock_file_path)

        def load(path_args: tuple[Any, ...]) -> None:
            file_path = get_path(path_args)
            if os.path.exists(file_path):
                with open(file_path, "rb") as cache_file:
                    loaded_cache[path_args] = pickle.load(cache_file)
            else:
                loaded_cache[path_args] = {}

        atexit.register(save)
        return persistent_cache_wrapper

    return persistent_cache_decorator


def update(cache_path_a: str, cache_path_b: str) -> None:
    """Cache a gets updated with cache b.

    In conflicts, values from b persist.
    B remains unchanged.
    """
    with open(cache_path_a, "rb") as cache_a:
        cache = pickle.load(cache_a)
    with open(cache_path_b, "rb") as cache_b:
        cache.update(pickle.load(cache_b))
    with open(cache_path_a, "wb") as cache_a:
        pickle.dump(cache, cache_a)


@persistent_cache((1, 0))
def test(tx: str, num: int = 1):
    print("Running function")
    return tx * num
