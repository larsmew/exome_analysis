# Conda for Bioinformatic Pipelines

## Table of contents

<!--ts-->
- [Install Conda package manager](#install-conda-package-manager)
- [Restore](#restore)
- [Restore](#Restore)
- [Install Conda package manager](Install Conda package manager)
<!--te-->

This document describes some simple commands to get the conda package up and running on 64-bit Linux machines (and likely MacOS machines).

## Install Conda package manager

There exist two versions of the conda package manager: **Miniconda** and **Anaconda**.

The difference is that Anaconda includes 150+ open source packages. If you are mostly doing bioinformatic pipelines I suggest you install **Miniconda** as we have to install most tools anyway.

### Install Miniconda

Run the following two commands to install conda. Follow the on-screen installation instructions.
The default location for conda is `~/miniconda3`.

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Remember to prepend the miniconda path to PATHS. Either, during install or by running the following:

Bash: `export PATH=~/miniconda3/bin:$PATH`

Fish: `set -U fish_user_paths ~/miniconda3/bin $fish_user_paths`

**Restart the shell for the changes to take effect.**

### Install Anaconda

To be done...

## Search for packages to install

Conda searches the anaconda repository as default

```bash
conda search <package_name>
```

**Example:**
```bash
conda search python
```

### Search BioConda

Most of the tools we are interested in will be located in the bioconda repository. To search this, we use:
```bash
conda search -c bioconda <package_name>
```

**Example:**
```bash
conda search -c bioconda samtools
```

## Install packages

In general:

```bash
conda install <package_name>
```

**Example:**
```bash
conda install python
```

### Install from BioConda

Installing from the bioconda repository is as easy as searching it:
```bash
conda install -c bioconda <package_name>
```

**Example:**
```bash
conda install -c bioconda samtools
conda install -c bioconda bwa
conda install -c bioconda picard
conda install -c bioconda gatk
conda install -c bioconda snakemake
```

## Backup and Restore conda environments

For reprodibility it can sometimes be necessary rerun analysis. Therefore it can be crucial to run the analysis exactly as it was done the first time with same set of packages and versions.

### Backup

To backup the current conda environment to a file use:

```bash
conda list --export > <file_name.txt>
```

**Example:**
```bash
conda list --export > conda_packages.txt
```

### Restore

To restore a conda environment without touching the current installation use the following command which creates a new conda environment with the specific packages defined from the backup file:

```bash
conda create --name <env> --file <file_name.txt>
```

Now, this environment can be activated:

```bash
source activate <env>
```

or deactivated:

```bash
source deactivate
```

**Example:**

```bash
conda create --name old_env --file conda_packages.txt
source activate old_env
```

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Cras ex risus, venenatis eget rhoncus eu, ullamcorper nec risus. Suspendisse auctor eros tortor. Pellentesque mattis cursus eros. Nulla id fermentum purus. Donec vulputate libero gravida urna egestas, non efficitur enim convallis. Vivamus accumsan magna vitae massa iaculis porta. Sed fringilla risus ac leo posuere condimentum. Ut commodo, nisl eu vehicula pellentesque, libero urna luctus ipsum, egestas lobortis neque tortor nec turpis. Phasellus est nulla, pulvinar id ullamcorper vel, efficitur ut mi. Curabitur elementum et nibh et fermentum.

Phasellus tincidunt pellentesque nibh, id blandit urna pellentesque ut. Maecenas hendrerit venenatis vestibulum. Integer at mauris tincidunt, mollis dui nec, sollicitudin risus. Mauris sed lectus et ex lobortis vestibulum. Morbi iaculis metus vel ipsum condimentum, at lobortis enim malesuada. Ut eu ex quis nisi ultrices posuere. Vivamus eu ultrices nisl. Proin sit amet eros accumsan nulla interdum porta in in leo. Vivamus sit amet urna diam. Nulla a enim ultricies, pharetra sapien id, sollicitudin tortor. Cras luctus varius nunc ut ullamcorper. Nullam sollicitudin rhoncus ante, sed hendrerit tellus lobortis ac.

Suspendisse et laoreet dolor. Praesent a tortor tincidunt, varius quam non, feugiat massa. Integer sollicitudin facilisis velit et gravida. Mauris ut dapibus enim. Cras felis erat, venenatis sit amet augue in, maximus tempus nisi. In et cursus diam, in mollis mi. Maecenas ut erat non ligula maximus scelerisque. Cras iaculis ornare finibus. Nullam gravida tellus nec augue placerat interdum. Aenean non pulvinar dui.

Pellentesque tincidunt euismod ultrices. Duis semper sit amet elit laoreet tristique. Donec tincidunt velit quis est sodales congue. Cras auctor scelerisque ex. Aenean luctus hendrerit sodales. Vivamus nec massa purus. Praesent pretium, ligula posuere laoreet scelerisque, tellus lorem aliquet dolor, vel ullamcorper mauris libero commodo augue. Ut purus odio, consequat eget ornare quis, euismod quis ligula. Phasellus nec mattis dui, laoreet efficitur ipsum. Morbi auctor ullamcorper ante a tempus. Maecenas viverra, odio sed tincidunt placerat, eros erat blandit lacus, in semper lectus nulla eu dolor. Mauris eu tortor purus.

Duis in nisi ut lorem sagittis volutpat. Quisque cursus augue sit amet turpis commodo commodo. Etiam ut blandit elit. Morbi ac est finibus, vestibulum ligula eget, imperdiet ipsum. Etiam lobortis leo vitae eros posuere blandit. Praesent mollis aliquam libero in vehicula. Integer ut lectus luctus, gravida elit dictum, facilisis elit. Vestibulum ut lectus nisl. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque nisl ex, posuere placerat justo imperdiet, molestie euismod neque. Maecenas ultricies quam in finibus condimentum.

Pellentesque posuere metus a mi efficitur, quis semper risus posuere. Etiam eget vulputate tortor. Proin sit amet ligula accumsan, tincidunt nisi eget, rutrum odio. Fusce tempus justo sed pharetra pretium. In nisl nisl, consequat ut dui a, rutrum blandit mauris. Nullam blandit mauris ex, nec malesuada nunc malesuada eget. Sed pharetra tristique commodo. Nulla consequat dignissim quam in feugiat. Duis imperdiet blandit diam, sit amet viverra nisl. Duis in quam at ligula accumsan rutrum vitae ac enim. Fusce fringilla sapien a purus hendrerit dapibus. Sed eget mollis orci.

Sed lacinia mollis maximus. Suspendisse pretium fermentum ante sed vehicula. Donec lobortis mi at libero venenatis vehicula. Quisque tempus metus justo, at fermentum erat dignissim vitae. In erat augue, elementum a leo quis, rutrum sagittis velit. Sed sit amet lacinia ligula. Mauris eget tempor nisi, a condimentum felis. Praesent tincidunt ex ut mi condimentum, volutpat porttitor est luctus. Duis tincidunt tincidunt consectetur. Nullam et ante quis diam consequat ornare. Curabitur ullamcorper mauris ac ex placerat, at vehicula orci lacinia. Nam ut metus sit amet lorem iaculis lacinia. Integer aliquet non eros sed tincidunt. Sed maximus orci nec lacus aliquet pulvinar. Etiam cursus lobortis ipsum, vel porttitor augue iaculis vel. Vivamus vitae tortor a nisi sagittis tempor ac et mauris.

Duis imperdiet gravida nibh. Aliquam erat volutpat. Cras posuere velit quis massa fermentum tincidunt. Donec at convallis mauris. Donec eget purus elementum, interdum metus nec, suscipit arcu. Morbi consequat iaculis ligula vel pellentesque. Nulla at orci sed magna efficitur condimentum in in massa. Suspendisse auctor orci quam, sed pulvinar mauris scelerisque eget. Sed laoreet massa at dolor lacinia porta. Vestibulum condimentum metus in elit imperdiet, eget pellentesque lorem posuere. Nulla sed neque pellentesque, aliquet elit et, pharetra tortor. Maecenas quis nunc in velit pulvinar viverra nec sit amet est.

Nunc auctor et neque sit amet vehicula. Sed enim velit, lacinia quis consequat ac, tempus in purus. Integer vulputate sapien quis lorem congue, sit amet vestibulum enim finibus. Etiam mi ipsum, fringilla nec iaculis sed, mattis a mauris. Maecenas cursus est non lectus molestie, at elementum urna aliquam. Nunc suscipit gravida est, sit amet blandit dolor euismod sed. In venenatis porta sapien. Phasellus odio leo, pharetra sit amet fringilla gravida, sagittis non dolor. Nunc sagittis mattis ultrices. Quisque eu metus eget nibh consectetur ultricies vitae ac nunc. Curabitur placerat dolor in erat vestibulum, laoreet accumsan tortor tempus. Curabitur varius convallis consectetur. Sed porta purus et libero sagittis condimentum.

Duis elit dolor, rutrum sit amet bibendum sit amet, rhoncus nec lectus. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Ut ligula odio, interdum fermentum odio at, commodo luctus nunc. Nulla aliquet lobortis pretium. Morbi luctus eros sed nunc tincidunt mattis. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Vivamus maximus pretium faucibus. Aliquam egestas quis lacus sodales blandit. Nullam dignissim nisl tellus, eu placerat mi imperdiet at. Maecenas euismod velit non nisl maximus venenatis. Sed vel ligula sodales dui auctor eleifend et sed leo. Mauris sed velit placerat, hendrerit ante pellentesque, ullamcorper lacus. Mauris id odio non metus varius posuere posuere a ex. Proin maximus, ligula in condimentum sodales, enim nunc ultricies tellus, a ornare sapien leo ac velit.

Etiam ornare tellus et varius sollicitudin. Duis gravida lorem eu lectus gravida, iaculis ornare massa pellentesque. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nam suscipit orci ac aliquam euismod. Ut sit amet elementum orci. Nunc interdum convallis blandit. Cras suscipit augue enim, eu consectetur dolor pellentesque in. Donec pretium faucibus diam, eu efficitur neque volutpat id. Sed vestibulum ut augue vel ultrices. Fusce dapibus lacus ut elit dapibus venenatis.

Vivamus nec diam luctus, rutrum ipsum non, placerat elit. Vestibulum tincidunt lectus enim, sit amet varius dolor convallis vel. Maecenas libero mi, dignissim ut nunc in, volutpat bibendum erat. Vivamus non euismod ipsum. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia curae; Aenean non consequat nisl, sed finibus mi. Duis blandit porttitor faucibus. Nam vulputate interdum malesuada. Nunc congue ornare massa et vulputate. Praesent finibus arcu at ipsum eleifend congue. Nunc ac elit purus. Aenean porta massa nec tellus tincidunt fermentum. Donec sed consequat urna. Vestibulum arcu erat, condimentum auctor maximus et, pellentesque non libero. Vivamus bibendum vulputate orci, vitae eleifend dui ullamcorper nec. Sed eu augue ut tellus ultricies porta.

Curabitur elementum quam quis felis pulvinar tristique. Nunc tincidunt orci eu risus rutrum euismod. Nullam neque tellus, mattis a leo at, consectetur volutpat elit. Praesent ullamcorper, nisi vel dictum condimentum, lectus nulla ornare orci, sit amet lobortis metus nunc quis elit. Fusce id massa non dui eleifend porta posuere vel libero. Fusce eget elit eu nisl consectetur tempor quis at nibh. Donec et imperdiet odio, non convallis nunc. Proin dapibus ipsum arcu, sed dictum magna maximus eget. Vivamus ultricies urna at odio rutrum fringilla. Pellentesque eget velit sapien. Sed ut fringilla tortor. Vestibulum cursus et quam eu imperdiet.

Nunc eu nisl quis ipsum pretium aliquam. Etiam leo dolor, aliquam nec ex ut, facilisis tristique tortor. Nunc at egestas augue. Donec convallis risus sed nunc vulputate fringilla. Donec dapibus in odio id porta. In hac habitasse platea dictumst. Proin posuere sem justo, nec maximus turpis ultricies eleifend. Pellentesque ac auctor velit, nec feugiat dui. Duis aliquet nibh id mauris auctor, nec vestibulum orci ultricies. Praesent enim ipsum, tempor fringilla odio sit amet, mattis ultricies libero. Etiam eget magna convallis nibh imperdiet dictum.

Fusce a faucibus dolor. Duis tristique interdum diam, ut egestas diam tempor nec. Morbi vestibulum porta mi eget fringilla. Aenean a est sem. Donec vel varius felis, ut placerat nisl. Vivamus sollicitudin egestas eleifend. Aenean commodo augue nisl, id euismod arcu vulputate et. Quisque in volutpat magna. Aenean et sodales ante, non porttitor ante.

Pellentesque convallis velit vel efficitur finibus. Lorem ipsum dolor sit amet, consectetur adipiscing elit. In sagittis purus quis lobortis vulputate. Vivamus vulputate rutrum orci. Etiam maximus elementum molestie. Nunc nisi elit, bibendum vitae auctor ut, tincidunt in purus. Sed et sodales ipsum, sed sodales est. Sed molestie odio quis ante bibendum egestas. Donec egestas eleifend erat, et faucibus est ultrices et. Nullam in tincidunt tellus. Donec ac elit ante. Cras ut sem massa. Fusce facilisis est eget ultricies maximus. Nunc congue eu ex a tempus.

Vivamus vulputate est ut congue rhoncus. Donec ac elementum ligula, at malesuada eros. Aenean et leo sit amet est elementum tincidunt a a lectus. Vestibulum a nibh neque. Curabitur posuere risus ut sodales accumsan. Etiam aliquam, erat id condimentum elementum, lorem enim feugiat ante, vitae maximus augue tellus a dolor. Quisque vel lectus purus. Suspendisse sit amet velit nibh. Morbi a diam id arcu eleifend lobortis. Aenean in porta lorem. Aliquam gravida scelerisque lorem, sit amet convallis magna auctor nec.

Maecenas leo lectus, bibendum non vehicula nec, fermentum a lorem. Phasellus tempus lacinia metus non efficitur. Pellentesque rutrum turpis rutrum tortor rhoncus hendrerit. Ut eget sapien elit. Vestibulum lacinia lorem ac imperdiet consequat. Donec vestibulum molestie sem at accumsan. Sed mattis faucibus sem quis tempus.

Fusce posuere fringilla sem, at elementum nibh. Pellentesque venenatis erat vitae leo accumsan, nec rhoncus ipsum maximus. Ut cursus lacus eu dui ultricies tempus. Suspendisse libero mauris, rhoncus sed sapien ut, egestas rutrum sem. Vivamus mi lacus, suscipit quis ipsum in, rhoncus blandit odio. Curabitur sapien metus, faucibus a sodales at, finibus ut sapien. Aenean vitae neque sed nisi facilisis lacinia eu quis magna. Maecenas semper pulvinar vulputate. Nam sed egestas sapien, vitae varius lectus. Phasellus finibus eu magna non accumsan. Morbi et magna nec odio interdum sodales. Integer pulvinar felis congue mi cursus, at vulputate odio euismod. Vestibulum et mattis elit. Nam posuere mollis ex, eget tincidunt odio sagittis ut. Aenean elementum sem porta tincidunt aliquam.

Pellentesque habitant morbi tristique senectus et netus et malesuada fames ac turpis egestas. Morbi porta neque ac vulputate viverra. Nunc arcu dolor, condimentum sit amet euismod quis, ultrices non risus. Quisque elit libero, tincidunt ut suscipit vitae, gravida at libero. Aenean sed ipsum lorem. Ut massa risus, blandit condimentum dapibus et, finibus ac nulla. Praesent porta egestas aliquam. Pellentesque habitant morbi tristique senectus et netus et malesuada fames ac turpis egestas. Suspendisse sed faucibus justo, nec pulvinar orci. In hac habitasse platea dictumst. Suspendisse auctor ac metus nec commodo. Mauris sed orci sit amet lacus sagittis iaculis at non sem. In hac habitasse platea dictumst. Nullam pulvinar cursus dui, ac scelerisque mi dictum nec.



## Restore
hej med dig 