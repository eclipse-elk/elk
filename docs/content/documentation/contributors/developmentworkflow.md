---
title: "Development Workflow"
menu:
  main:
    identifier: "DevelopmentWorkflow"
    parent: "ELKContributors"
    weight: 20
---

Once you [have everything set up]({{< relref "documentation/contributors/developmentsetup.md" >}}) it is time to start getting your hands at the code. The exact steps depend on how you have set up your development environment and your taste in things (and stuff), but here's what you will generally want to do.

1. **Set Up an Eclipse Installation**

    This used to be rather hard, but thankfully people have seen that problem and started to do something about it. The result is the _Oomph Eclipse Installer_. The easiest way now is to go [to the Eclipse Downloads page](https://www.eclipse.org/downloads/index.php) and click the prominent _Download_ button there. Execute the installer and [follow this guide]({{< relref "documentation/contributors/developmentworkflow/installingwithoomph.md" >}}) (includes pictures!).

1. **Find a Git Repository to Work With**

    While working you will probably want to create a bunch of commits and push them to a repository. However, our installer by default simply clones our official ELK repository at GitHub, which you don't have permission to write to. What you will want to do, thus, is to get your hands on another GitHub repository that you do have permission to write to. The best way is to simply fork our main repository by clicking the _Fork_ link at the upper right corner of [our repository page](https://github.com/eclipse-elk/elk). Add your new repository as a second remote to your local clone and push your commits there.

1. **Tickets**

    Of course we have [a list of issues](https://github.com/eclipse-elk/elk/issues). If you want to work on one of them, you should mention your intentions there to let us know that we should not start working on the issue ourselves. Also, we might have helpful input for you to think about. If you want to implement a feature that we do not have an issue for yet, open one.

1. **Proper Commit Messages**

    This is an example of what we consider to be a proper commit message:

    ```plain
    Core, Layered: Fix labels with layout directions #58

    Using node labels inside compound nodes that are configured to use
    different layout directions resulted in problems with calculated insets.

    Signed-off-by: Christoph Daniel Schulze <cds@informatik.uni-kiel.de>
    ```

    It consists of three parts:

    1. A first line that summarizes the commit. The summary starts by mentioning the general parts of the project that the commit affects, followed by a summary and a reference to the issue the commit is associated with. If there is no such issue, you can of course omit that.

    1. An optional paragraph (separated from the first line by a blank line) that provides more information, if necessary.

    1. A sign-off, preceded by a blank line. This is your name followed by the e-mail address that belongs to your [Eclipse.org](https://www.eclipse.org/) account. The sign-off is important. With this, you confirm that to the best of your knowledge your commit complies with the [Certificate of Origin](http://www.eclipse.org/legal/CoO.php) that you promised to honor when you signed the Eclipse ECA.

1. **Contributing Your Code**

    To actually contribute your code to the ELK project, go [to our pull requests page](https://github.com/eclipse-elk/elk/pulls) and click _New pull request_. Since you want us to pull in changes from your fork, click the _Compare across forks_ link and select the branch that contains the relevant changes. We will then review your code. Chances are that we will reject your pull request on your first try and ask you to make some changes. Don't let that discourage you (also, we try to be nice about it). Simply make the changes and create another pull request.
