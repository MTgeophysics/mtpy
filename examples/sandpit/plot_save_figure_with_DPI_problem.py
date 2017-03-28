"""
Description:
save plot behaves differently when dpi= is specified.
https://github.com/matplotlib/matplotlib/issues/7820/
Author: fei.zhang@ga.gov.au

Date: 2017-03-15
"""

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns

def demo_example():
    """
    this need http access of data
    :return:
    """
    df = pd.read_csv('http://www.who.int/childgrowth/standards/wfa_boys_p_exp.txt', sep='\t', index_col=0)

    DAYS_IN_MONTH = 30.4375
    DAYS_IN_YEAR = DAYS_IN_MONTH * 12

    percentile_list = []

    for i in range(7):
        percentile_list.append((tuple(df.columns[[3 + i,17 - i]])))

    df.index = df.reset_index()['Age'].div(DAYS_IN_MONTH)
    df = df.loc[df.index < 12]

    fig, ax = plt.subplots()
    df['P50'].plot(ax=ax, color='w')

    x = df.index

    color_list = sns.color_palette("Blues", len(percentile_list))

    for i, j in enumerate(percentile_list):
        y1 = df[j[0]]
        y2 = df[j[1]]
        ax.fill_between(x, y1, y2, color=color_list[i])

    text_annotations= []

    for i in df.columns[3:]:
        annotation = ax.annotate(i, xy=(12, df.loc[:,i].iloc[-1]),
                                 xytext=(5,0), textcoords="offset points",
                                 va="center", ha="left", size=8)
        text_annotations.append(annotation)

    fig.savefig('weight_chart.png', bbox_inches='tight')
    fig.savefig('weight_chart_dpi.png', bbox_inches='tight', dpi=300)


def a_simpler_example():

    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [1, 2, 3])

    ax.annotate('annotation 1', xy=(3, 1.5))
    ax.annotate('annotation 2', xy=(3, 2.5))

    fig.savefig('chart_dpiNone.png', bbox_inches='tight')
    fig.savefig('chart_dpi100.png', bbox_inches='tight', dpi=200)
    fig.savefig('chart_dpi200.png', bbox_inches='tight', dpi=200)
    fig.savefig('chart_dpi300.png', bbox_inches='tight', dpi=300)

if __name__=="__main__":
    a_simpler_example()
