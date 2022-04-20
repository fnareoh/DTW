import matplotlib.pyplot as plt
import random
import time

sizes = [i * 10 ** 5 for i in range(50)]
runtimes = []
for size in sizes:
    s = set(range(size))
    t = set(range(0, size, 2))
    # Start track time ...
    t1 = time.time()
    s.difference(t)
    t2 = time.time()
    # ... end track time

    runtimes.append(t2 - t1)
plt.plot(sizes, runtimes)
plt.ylabel("Runtime (s)")
plt.xlabel("Set Size")
plt.show()
