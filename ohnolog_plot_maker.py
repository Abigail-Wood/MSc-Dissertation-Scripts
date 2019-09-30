import matplotlib
import matplotlib.pyplot as plt
import numpy as np


labels = ['Innate hits', 'InnateDB', 'Innate excluded', 'All hits', 'All human genes']
singletons = [11.9,11.2,9.77,23.0,7.38]
ssds = [47.1,58.9,63.4,46.8,67.5]
ohnologs = [40.9,29.9,26.9,30.2,25.1]

x = np.arange(len(labels))  # the label locations
width = 0.6  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/3, singletons, width/3, label='Singleton')
rects2 = ax.bar(x, ssds, width/3, label='SSD')
rects3 = ax.bar(x + width/3, ohnologs, width/3, label='Ohnolog')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Paralog status (% of dataset, to 3 s.f.)')
ax.set_title('Paralog status')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 3, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)
autolabel(rects3)

fig.tight_layout()

plt.show()
