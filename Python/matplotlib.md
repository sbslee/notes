# matplotlib

## Table of contents

* [Frequently used commands](#Frequently-used-commands)
* [Combining subplots](#Combining-subplots)
* [Setting space between subplots](#Setting-space-between-subplots)
* [Setting figure style globally](#Setting-figure-style-globally)
* [Setting figure style temporarily](#Setting-figure-style-temporarily)

## Frequently used commands <a name="Frequently-used-commands"></a>

To set figure title:

```
fig.suptitle('This is a somewhat long figure title')
```

To set figure title in tight layout:

```
fig.suptitle('This is a somewhat long figure title')
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
```

To remove a subplot:

```
ax.clear()
ax.axis('off')
ax.set_visible(False)
ax.remove()
```

To set widths and heights of suplots:

```
fig, [ax1, ax2] = plt.subplots(1, 2, gridspec_kw={'width_ratios': [9, 1], 'height_ratios': [1, 3]})
```

To remove legend title:

```
plt.gca().legend().set_title('')
```

To set default figure style:

```
matplotlib.rc_file_defaults()
```

To remove gaps between subplots:

```
plt.subplots_adjust(wspace=0, hspace=0)
```

To set font sizes:

```
ax.xaxis.label.set_size(20)
ax.yaxis.label.set_size(20)
ax.tick_params(axis='x', which='major', labelsize=15)
ax.tick_params(axis='y', which='major', labelsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_title('My subplot title', fontsize=30)
ax.legend(handles, labels, fontsize=20, title=title, title_fontsize=25)
```

## Combining subplots <a name="Combining-subplots"></a>

```
# Source: https://matplotlib.org/3.1.1/gallery/subplots_axes_and_figures/gridspec_and_subplots.html

import matplotlib.pyplot as plt

fig, axs = plt.subplots(ncols=3, nrows=3)
gs = axs[1, 2].get_gridspec()
# remove the underlying axes
for ax in axs[1:, -1]:
    ax.remove()
axbig = fig.add_subplot(gs[1:, -1])
axbig.annotate('Big Axes \nGridSpec[1:, -1]', (0.1, 0.5),
               xycoords='axes fraction', va='center')

fig.tight_layout()

plt.show()
```


## Setting space between subplots <a name="Setting-space-between-subplots"></a>

```
# Source: https://stackoverflow.com/questions/49781442/matlibplot-how-to-add-space-between-some-subplots

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# Simple data to display in various forms
x = np.linspace(0, 2 * np.pi, 400)
y = np.sin(x ** 2)

f = plt.figure(figsize=(10,10))
gs0 = gridspec.GridSpec(2, 1)

gs00 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[0], hspace=0)
ax0 = f.add_subplot(gs00[0])
ax0.plot(x, y)
ax0.set_title('Panel: A')
ax1 = f.add_subplot(gs00[1], sharex=ax0)
ax1.plot(x, y**2)

gs01 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[1], hspace=0)
ax2 = f.add_subplot(gs01[0])
ax2.plot(x, y**3)
ax2.set_title('Panel: B')
ax3 = f.add_subplot(gs01[1], sharex=ax0)
ax3.plot(x, y**4)

plt.show()
```

## Setting figure style globally <a name="Setting-figure-style-globally"></a>

```
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline

plt.style.use('ggplot')

plt.plot(np.sin(np.linspace(0, 2 * np.pi)), 'r-o')
```

If you want to use the `seaborn` package's default style:

```
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
%matplotlib inline

sns.set()

plt.plot(np.sin(np.linspace(0, 2 * np.pi)), 'r-o')
```

## Setting figure style temporarily <a name="Setting-figure-style-temporarily"></a>

```
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline

with plt.style.context('ggplot'):
    plt.plot(np.sin(np.linspace(0, 2 * np.pi)), 'r-o')
```

If you want to use the `seaborn` package's default style:

```
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
%matplotlib inline

with sns.axes_style('darkgrid'):
    plt.plot(np.sin(np.linspace(0, 2 * np.pi)), 'r-o')
```