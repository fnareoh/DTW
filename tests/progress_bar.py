import sys


def update_progress(progress):
    barLength = 50  # Modify this to change the length of the progress bar
    status = ""
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength * progress))
    text = "\rPercent: [{0}] {1}% {2}".format(
        "#" * block + "-" * (barLength - block), round(progress * 100, 2), status
    )
    sys.stderr.write(text)
    sys.stderr.flush()


if __name__ == "__main__":
    n = 10000000
    for i in range(n):
        if i % 1000 == 0:  # Avoid to update at each step (time consuming)
            update_progress(i / float(n))
            # your real code here
    update_progress(1)  # Ends the line & writes "Done"
