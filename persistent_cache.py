# Implements a decorator, that caches function call values cross runtimes

from typing import Any
import atexit
import pickle
import os.path
import os
from time import sleep


print("Using persistent caching provided by persistent_cache.py")


def retry(function):
    """Retries a function until no exception occurs within it's execution.

    Stops after 10 000 attempts with 0.1 seconds of delay in between.
    In this case the function is called a final time and any exception is
    carried on.
    """
    timeout = 0.1
    retries = 10_000

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
        print(
            f"Function {function.__name__} could not be executed in {retries} attempts."
        )
        return function(*args)


def persistent_cache(
    path_arg_mask: tuple[int, ...] = (),
    filename: str | None = None,
):
    """A cache that saves it's contents to disks and can reuse them on the next run.

    IMPORTANT: All arguments of the decorated function must be hashable.

    The cache is saved under the name of the function. Multiple functions of
    the same name may save into the same file/folder. To prevent this, identify
    provide a unique enough identifier, e.g. the filename of the python script
    that defines the function.

    Args:
        path_arg_mask: Tuple with integers in the range 0 to #function arguments -1.
            The argument of the function with an index in this tuple will be used as
            a path argument.
            Function values with different path arguments will
            be saved in different files. The path to the save-file is named after
            these path arguments (in the order in which they appear in the tuple).
            Therefore arguments whose string representation is not a valid path
            may result in exceptions.
        filename: A string, which defines the first part (first folder) of the
            save path. Usefull for scoping of the caches since they are only
            identified by the function name by default. If unique enough it is
            recommended to use the filename of the python script that defines the
            function.
    """
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
        def lock(file_path: str):
            """Waits, until a file is no longer locked and locks it.

            Creates a lock file (<filepath>.lock).
            If it already exists, waits for it to be deleted (i.e. the file
            being unlocked) before locking the file again.

            Note that this system only works by convention. Other programms or a user
            may change the file while being locked.

            Note also that if the file takes to long to unlock (>1000 seconds), this
            will end in an error thrown.

            Args:
                file_path: The path of the file to be locked.
            """
            with open(file_path + ".lock", "x"):
                # Using open with the x argument creates a file but raises an exception
                # if the file already exists.
                pass

        def unlock(file_path: str):
            """Unlocks the file again.

            I.e. deletes the lock file.

            Args:
                file_path: the path of the locked file.
            """
            os.remove(file_path + ".lock")

        def save():
            # save files
            # gets called when the programm is terminated.

            for path_args, sub_cache in loaded_cache.items():
                # create dir if it does not exists
                file_path = get_path(path_args)
                dir_path = os.path.dirname(file_path)
                if not os.path.isdir(dir_path):
                    os.makedirs(dir_path)
                # lock file
                lock(file_path)
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
                unlock(file_path)

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

    Args:
        chache_path_a: Path to cache a
        chache_path_b: Path to cache b
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
