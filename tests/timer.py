import time


# This class allows us to monitor the time spent by a function
class Timer:
    # Function called right before the execution of the function
    def __enter__(self):
        self.t1 = time.perf_counter()
        return self

    # Function called right after the function
    def __exit__(self, type, value, traceback):
        self.t2 = time.perf_counter()
        self.t = self.t2 - self.t1

    # Function that prints on the shell the time spent by the instructions
    def print(self, template: str = "{}"):
        print(template.format(round(self.t, 2)))

    def value(self):
        return self.t


if __name__ == "__main__":
    # Exemple how to use this class
    with Timer() as total_time:  # time all instructions in the ’with’ statements
        for i in range(5):
            time.sleep(0.471)
    total_time.print("Durée = {} secondes")
