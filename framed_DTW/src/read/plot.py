import matplotlib.pyplot as plt
from statistics import mean, stdev
from sys import argv
from glob import glob
import csv


def reorder(x, y_dict):
    sorted_y = dict()
    for k, y in y_dict.items():
        zipped_lists = zip(x, y)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        _, sorted_y[k] = [list(tuple) for tuple in tuples]
    return sorted_y


def parse_csv(filename):
    output = dict()
    with open(filename, "r") as file:
        reader_csv = csv.reader(file)
        reader = [row for row in reader_csv]
        details = reader[0]
        assert details[4] == "proba insertion Deletion Substitution"
        output["IDS"] = float(details[5])
        assert details[6] == "proba Homopolymer extension"
        output["H"] = float(details[7])
        output["IDS_distrib"] = details[8]
        assert reader[1][-1] == "min ED" and reader[1][-2] == "min DTW"
        output["ED"] = [int(row[-1]) for row in reader[2:]]
        output["DTW"] = [int(row[-2]) for row in reader[2:]]
        output["mean_ED"], output["stdev_ED"] = mean(output["ED"]), stdev(output["ED"])
        output["mean_DTW"], output["stdev_DTW"] = mean(output["DTW"]), stdev(
            output["DTW"]
        )
        assert reader[1][1] == "nb_IDS" and reader[1][2] == "nb_H"
        output["nb_IDS"] = [int(row[1]) for row in reader[2:]]
        output["mean_IDS"], output["stdev_IDS"] = mean(output["nb_IDS"]), stdev(
            output["nb_IDS"]
        )
        output["nb_H"] = [int(row[2]) for row in reader[2:]]
    return output


def parse_args():
    other = {"IDS": "H", "H": "IDS"}
    param = {"IDS": "*", "H": "*"}
    param["basename"] = argv[1]
    tmp_parse = param["basename"].split("_")
    assert tmp_parse[-2] == "N"
    param["N"] = int(tmp_parse[-1])
    param["fixed"] = argv[2]
    assert param["fixed"] in ["IDS", "H"]
    param["varying"] = other[param["fixed"]]
    param[param["fixed"]] = float(argv[3])
    param["IDS_distrib"] = None

    # Get values
    x = []
    y_mean = {"DTW": [], "ED": [], "IDS": []}
    y_stdev = {"DTW": [], "ED": [], "IDS": []}
    all_values = {"IDS": [], "ED": [], "DTW": []}
    for name in glob(
        f"{param['basename']}_IDS_{param['IDS']}_H_{param['H']}_IDS_distrib_*.csv"
    ):
        output = parse_csv(name)

        if (
            param["IDS_distrib"] is None
            or param["IDS_distrib"] == output["IDS_distrib"]
        ):
            param["IDS_distrib"] = output["IDS_distrib"]
            x.append(output[param["varying"]])
            all_values["IDS"] += output["nb_IDS"]
            all_values["ED"] += output["ED"]
            all_values["DTW"] += output["DTW"]
            y_mean["DTW"].append(output["mean_DTW"])
            y_stdev["DTW"].append(output["stdev_DTW"])
            y_mean["ED"].append(output["mean_ED"])
            y_stdev["ED"].append(output["stdev_ED"])
            y_mean["IDS"].append(output["mean_IDS"])
            y_stdev["IDS"].append(output["stdev_IDS"])

    x, y_mean, y_stdev = sorted(x), reorder(x, y_mean), reorder(x, y_stdev)
    plot_avg_mean_stdev(param, x, y_mean, y_stdev)
    plt.figure().clear()

    x = all_values["IDS"]
    x, y = sorted(x), reorder(x, all_values)
    simple_plot(param, x, y)


def simple_plot(param, x, y):
    x_legend = {
        "IDS": f"Number of  IN, DEL or SUB error with repartition {param['IDS_distrib']}",
    }
    for k in y:
        plt.plot(x, y[k], label=k)
    plt.ylabel("Distance")
    plt.xlabel(x_legend[param["varying"]])
    plt.title(
        f"Edit and dynamic time warping distance as a function\n of the number of insertion deletion and substitution"
    )
    plt.legend(loc="upper left")
    plt.savefig(
        f"{param['basename']}_nb_IDS_fixed_{param['fixed']}_{param[param['fixed']]}.png"
    )


def plot_avg_mean_stdev(param, x, y_mean, y_stdev):
    x_legend = {
        "IDS": "Probability of IN, DEL or SUB error",
        "H": "Probability of homopolymer extension error",
    }
    for k in y_mean:
        plt.errorbar(x, y_mean[k], y_stdev[k], label=k, fmt="-o")
    plt.ylabel("Average distance")
    plt.xlabel(x_legend[param["varying"]])
    plt.title(
        f"Average distance for fixed {x_legend[param['fixed']]}\n equal to {param[param['fixed']]}, each point averaged over {param['N']} sequences"
    )
    plt.legend(loc="upper left")
    plt.savefig(
        f"{param['basename']}_fixed_{param['fixed']}_{param[param['fixed']]}_IDS_distrib_{param['IDS_distrib']}.png"
    )


def main():
    parse_args()


if __name__ == "__main__":
    main()
