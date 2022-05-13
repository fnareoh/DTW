import matplotlib.pyplot as plt
from statistics import mean,stdev
from sys import argv
from glob import glob
import csv

def reorder(x, y_dict,y_stdev):
    sorted_y = dict()
    sorted_y_stdev = dict()
    for k,y in y_dict.items():
        zipped_lists = zip(x, y, y_stdev[k])
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        _,sorted_y[k],sorted_y_stdev[k] = [list(tuple) for tuple in tuples]
    return sorted(x),sorted_y,sorted_y_stdev

def parse_csv(filename):
    output = dict()
    with open(filename, 'r') as file:
        reader_csv = csv.reader(file)
        reader = [row for row in reader_csv]
        details = reader[0]
        assert details[4]=='proba insertion Deletion Substitution'
        output['IDS'] = float(details[5])
        assert details[6]=='proba Homopolymer extension'
        output['H'] = float(details[7])
        assert reader[1][-1]== 'min ED' and reader[1][-2]== 'min DTW'
        output['ED'] = [int(row[-1]) for row in reader[2:]]
        output['DTW'] = [int(row[-2]) for row in reader[2:]]
        output['mean_ED'],output['stdev_ED'] = mean(output['ED']), stdev(output['ED'])
        output['mean_DTW'],output['stdev_DTW'] = mean(output['DTW']), stdev(output['DTW'])
        assert reader[1][1]== 'nb_IDS' and reader[1][2]== 'nb_H'
        output['nb_IDS'] = [int(row[1]) for row in reader[2:]]
        output['nb_H'] = [int(row[2]) for row in reader[2:]]
    return output

def parse_args():
    other = {"IDS":"H","H":"IDS"}
    param={"IDS":"*","H":"*"}
    param['basename'] = argv[1]
    tmp_parse = param['basename'].split('_')
    assert tmp_parse[-2]=='N'
    param['N'] = int(tmp_parse[-1])
    param['fixed'] = argv[2]
    assert param['fixed'] in ["IDS","H"]
    param['varying'] = other[param['fixed']]
    param[param['fixed']] = float(argv[3])

    # Get values
    x = []
    y_mean = {'DTW':[],'ED':[]}
    y_stdev = {'DTW':[],'ED':[]}
    for name in glob(f"{param['basename']}_IDS_{param['IDS']}_H_{param['H']}.csv"):
        output=parse_csv(name)
        x.append(output[param['varying']])
        y_mean['DTW'].append(output['mean_DTW'])
        y_stdev['DTW'].append(output['stdev_DTW'])
        y_mean['ED'].append(output['mean_ED'])
        y_stdev['ED'].append(output['stdev_ED'])

    return param,*reorder(x,y_mean,y_stdev)

def plot(param,x,y_mean,y_stdev):
    x_legend = {'IDS':'Probability of IN, DEL or SUB error',
    'H':'Probability of homopolymer extension error'}
    for k in  y_mean:
        plt.errorbar(x, y_mean[k], y_stdev[k],label=k,fmt='-o')
    plt.ylabel('Average distance')
    plt.xlabel(x_legend[param['varying']])
    plt.title(f"Average distance for fixed {x_legend[param['fixed']]}\n equal to {param[param['fixed']]}, each point averaged over {param['N']} sequences")
    plt.legend(loc="upper left")
    plt.savefig(f"{param['basename']}_fixed_{param['fixed']}_{param[param['fixed']]}.png")

def main():
    param,x,y_mean,y_stdev = parse_args()
    plot(param,x,y_mean,y_stdev)

if __name__ == "__main__":
    main()


