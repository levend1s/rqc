import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import numpy as np

import pod5 as p5

# Using the example pod5 file provided
example_pod5 = "/Users/joshualevendis/Downloads/FAZ60852_d9c3cda9_6cfcc167_15.pod5"
selected_read_id = '945af20d-09eb-43b9-ba54-17d14ed43a29'

with p5.Reader(example_pod5) as reader:

    # Read the selected read from the pod5 file
    # next() is required here as Reader.reads() returns a Generator
    read = next(reader.reads(selection=[selected_read_id]))

    # Get the signal data and sample rate
    sample_rate = read.run_info.sample_rate
    signal = read.signal

    # Compute the time steps over the sampling period
    time = np.arange(len(signal)) / sample_rate

    print(time)
    print(signal)
    print(len(signal))

    # Plot using matplotlib
    # plt.plot(time, signal)
    # plt.show()
    N=len(signal)
    stop = time[-1]

    # figure preparation
    fig, ax = plt.subplots(1, 1, figsize = (8*0.9, 6*0.9))
    displayed_period = int(2*sample_rate)
    span = int(N/stop/sample_rate)

    def animation(i):
        # delete previous frame
        ax.cla()

        # plot and set axes limits
        ax.plot(time[span*i: 1 + span*(i + displayed_period)],
                signal[span*i: 1 + span*(i + displayed_period)])
        ax.set_xlim([time[span*i], time[span*(i + displayed_period)]])
        ax.set_ylim([1.1*signal.min(), 1.1*signal.max()])

    # run animation
    anim = FuncAnimation(fig, animation, frames = int(len(time)/span - 1), interval = 1)
    plt.show()