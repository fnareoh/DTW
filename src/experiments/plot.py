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
        assert details[4] == "percentage of Indel"
        output["ID"] = float(details[5])
        assert details[6] == "percentage of SNP"
        output["SNP"] = float(details[7])
        assert details[8] == "proba Homopolymer extension"
        output["H"] = float(details[9])
        assert details[10] == "proba substitution sequencing"
        output["seq_S"] = float(details[11])
        assert (
            reader[1][1] == "biological_var"
            and reader[1][2] == "seq_S"
            and reader[1][3] == "seq_H"
        )
        output["biological_var"] = [int(row[1]) for row in reader[2:]]
        output["mean_bio"], output["stdev_bio"] = mean(output["biological_var"]), stdev(
            output["biological_var"]
        )
        output["seq_S"] = [int(row[2]) for row in reader[2:]]
        output["nb_homopoly"] = [int(row[3]) for row in reader[2:]]
        assert reader[1][-1] == "min ED" and reader[1][-2] == "min DTW"
        output["ED"] = [int(row[-1]) for row in reader[2:]]
        output["DTW"] = [int(row[-2]) for row in reader[2:]]
    return output


def minus_bio_var(output):
    output["DTW_less_bio"] = [
        output["DTW"][i] - output["biological_var"][i]
        for i in range(len(output["DTW"]))
    ]
    output["ED_less_bio"] = [
        output["ED"][i] - output["biological_var"][i] for i in range(len(output["DTW"]))
    ]


def avg_same_x(x, y_dict):
    y_mean, y_stdev = {}, {}
    current_y = {}
    unique_x = []
    for k, y in y_dict.items():
        y_mean[k] = []
        y_stdev[k] = []
        current_y[k] = []
    last_x = 0
    for i in range(len(x)):
        if x[i] == last_x:
            for k, y in y_dict.items():
                current_y[k].append(y[i])
        else:
            unique_x.append(last_x)
            for k, y in y_dict.items():
                if len(current_y[k]) > 0:
                    y_mean[k].append(mean(current_y[k]))
                    if len(current_y[k]) > 1:
                        y_stdev[k].append(stdev(current_y[k]))
                    else:
                        y_stdev[k].append(0)
                current_y[k] = [y[i]]
            last_x = x[i]
    unique_x.append(last_x)
    for k, y in y_dict.items():
        if len(current_y[k]) > 0:
            y_mean[k].append(mean(current_y[k]))
            y_stdev[k].append(stdev(current_y[k]))
        current_y[k] = []
    return unique_x, y_mean, y_stdev


def parse_args():
    other = {"ID": "H", "H": "ID"}
    param = {"ID": "*", "H": "*", "SNP": 1}
    param["basename"] = argv[1]
    tmp_parse = param["basename"].split("_")
    assert tmp_parse[-2] == "N"
    param["N"] = int(tmp_parse[-1])
    param["fixed"] = argv[2]
    assert param["fixed"] in ["ID", "H"]
    param["varying"] = other[param["fixed"]]
    param[param["fixed"]] = round(float(argv[3]), 2)

    # Get values
    x = []
    y_mean = {"DTW": [], "ED": [], "bio": []}
    y_stdev = {"DTW": [], "ED": [], "bio": []}
    all_values = {
        "biological_var": [],
        "nb_homopoly": [],
        "ED": [],
        "DTW": [],
        "ED_less_bio": [],
        "DTW_less_bio": [],
    }
    for name in glob(
        f"{param['basename']}_ID_{param['ID']}_SNP_{param['SNP']}_H_{param['H']}_seqS_0.001*.csv"
    ):
        output = parse_csv(name)
        minus_bio_var(output)  # Deducts the bio distance !
        x.append(output[param["varying"]])
        all_values["biological_var"] += output["biological_var"]
        all_values["nb_homopoly"] += output["nb_homopoly"]
        all_values["ED"] += output["ED"]
        all_values["DTW"] += output["DTW"]
        all_values["ED_less_bio"] += output["ED_less_bio"]
        all_values["DTW_less_bio"] += output["DTW_less_bio"]
        y_mean["DTW"].append(mean(output["DTW_less_bio"]))
        y_stdev["DTW"].append(stdev(output["DTW_less_bio"]))
        y_mean["ED"].append(mean(output["ED_less_bio"]))
        y_stdev["ED"].append(stdev(output["ED_less_bio"]))
        y_mean["bio"].append(mean(output["biological_var"]))
        y_stdev["bio"].append(stdev(output["biological_var"]))

    x, y_mean, y_stdev = sorted(x), reorder(x, y_mean), reorder(x, y_stdev)
    print("mean DTW: ", y_mean["DTW"])
    print("stdev DTW: ", y_stdev["DTW"])
    plot_avg_mean_stdev(param, x, y_mean, y_stdev)
    plt.figure().clear()

    print(
        f"Average ratio DTW/bio_dist: {mean([all_values['DTW'][i]/(all_values['biological_var'][i]) for i in range(len(all_values['DTW'])) if all_values['biological_var'][i] > 0])}"
    )
    print(
        f"Stdev ratio DTW/bio_dist: {stdev([all_values['DTW'][i]/(all_values['biological_var'][i]) for i in range(len(all_values['DTW'])) if all_values['biological_var'][i] > 0])}"
    )
    print(
        f"Average ratio ED/bio_dist: {mean([all_values['ED'][i]/(all_values['biological_var'][i]) for i in range(len(all_values['DTW'])) if all_values['biological_var'][i] > 0])}"
    )
    print(
        f"Stdev ratio ED/bio_dist: {stdev([all_values['ED'][i]/(all_values['biological_var'][i]) for i in range(len(all_values['DTW'])) if all_values['biological_var'][i] > 0])}"
    )
    x = all_values["nb_homopoly"]
    all_values.pop("biological_var", None)
    all_values["DTW"] = all_values["DTW_less_bio"]
    all_values["ED"] = all_values["ED_less_bio"]
    all_values.pop("DTW_less_bio", None)
    all_values.pop("ED_less_bio", None)
    x, y = sorted(x), reorder(x, all_values)
    x, y_mean, y_stdev = avg_same_x(x, y)
    plot_avg_mean_stdev(param, x, y_mean, y_stdev, "_bio_var")


def simple_plot(param, x, y):
    x_legend = {
        "H": "Number of homopolymer extension",
        "bio": f"Number of  IN, DEL or SUB with repartition ",
    }
    for k in y:
        plt.plot(x, y[k], label=k)
    # plt.plot(x, [1] * len(x), color="black")
    plt.ylabel("Distance - biological variation")
    plt.xlabel(x_legend[param["varying"]])
    # plt.yscale("log")
    plt.title(
        f"Edit and dynamic time warping distance as a function\n of the number of homopolymer extensions"
    )
    plt.legend(loc="upper left")
    plt.savefig(
        f"{param['basename']}_nb_IDS_fixed_{param['fixed']}_{param[param['fixed']]}.png"
    )


def plot_avg_mean_stdev(param, x, y_mean, y_stdev, suffix=""):
    x_legend = {
        "ID": "Probability of IN, DEL or SUB error",
        "H": "Probability of homopolymer extension error",
    }
    for k in y_mean:
        # plt.plot(x, y_mean[k], label=k)
        if k != "bio":
            plt.errorbar(x, y_mean[k], y_stdev[k], label=k, fmt="-o")
    plt.ylabel("Average (distance - biological diversity)")
    plt.xlabel(x_legend[param["varying"]])
    # plt.title(
    #    f"Average distance for fixed {x_legend[param['fixed']]}\n equal to {param[param['fixed']]}, each point averaged over {param['N']} sequences"
    # )
    plt.legend(loc="upper left")
    plt.savefig(
        f"{param['basename']}_fixed_{param['fixed']}_{param[param['fixed']]}"
        + suffix
        + ".png"
    )


def main():
    parse_args()


if __name__ == "__main__":
    main()
