---
title: "Installing With Oomph"
menu:
  main:
    identifier: "InstallingWithOomph"
    parent: "DevelopmentWorkflow"
    weight: 10
---

Setting up an Eclipse to work on the _Eclipse Layout Kernel_ is comparatively easy when using the Oomph Installer. Follow this step-by-step guide for great glory.

1. Start the installer. This is what it will look like, more or less:

    {{< image src="oomph-install_step1.png" alt="Oomph Installation, Step 1" >}}

1. Click on the Hamburger menu (the three stripes at the top right). A menu with additional options will magically appear:

    {{< image src="oomph-install_step2.png" alt="Oomph Installation, Step 2" >}}

1. Click the _Advanced Mode_ button (we're programmers, we don't want to deal with the simple stuff). Here's what advanced mode looks like:

    {{< image src="oomph-install_step3.png" alt="Oomph Installation, Step 3" >}}

1. This page of the installer allows you to choose which basic Eclipse distribution you want to install. The _Eclipse Layout Kernel_ setup will use as the base upon which to add stuff you will need to hack ELK. Simply choosing the _Eclipse IDE for Java Developers_ and clicking _Next_ should work fine. You may want to make sure that you don't select _Latest Release_ in the _Product Version_ combo box at the bottom, but instead choose the latest release specifically. For example, instead of _Latest Release (Neon)_ select _Neon_ directly. This will keep Oomph from updating your installation and breaking things later on, which quite frankly will spare you the odd nervous breakdown.

    Anyway, after clicking the _Next_ button, this is what will appear:

    {{< image src="oomph-install_step4.png" alt="Oomph Installation, Step 4" >}}

1. It is now time to choose which project you actually want to work on. Find the _Eclipse Layout Kernel_ item (in the _Eclipse Projects_ category) and select it by clicking its check box to the left. Click _Next_ to see something like this:

    {{< image src="oomph-install_step5.png" alt="Oomph Installation, Step 5" >}}

1. Make sure everything looks good. You can usually leave the defaults, although you may need to change the _Eclipse Layout Kernel Github repository_ URL to whatever you see when you click on _Clone or download_ [on our repository page](https://github.com/eclipse-elk/elk). Clicking _Next_ will allow you to review what the installer wants to do. Unless you understand what that is, simply ignore it and click _Finish_.

1. Once the setup has finished the basic installation, your new Eclipse will open and proceed to clone the GitHub repository and to set everything up properly. Once all ELK projects appear in your workspace, you're set and ready to start developing.
