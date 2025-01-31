---
title: "0.8.1"
menu:
  main:
    parent: "ReleaseNotes"
    weight: -81
---

* [Release log](https://projects.eclipse.org/projects/modeling.elk/releases/0.8.1)
* [Documentation](https://download.eclipse.org/elk/updates/releases/0.8.1/elk-0.8.1-docs.zip)
* [Update site](https://download.eclipse.org/elk/updates/releases/0.8.1/)
* [Zipped update site](https://download.eclipse.org/elk/updates/releases/0.8.1/elk-0.8.1.zip) (for offline use)
* [Maven central](https://repo.maven.apache.org/maven2/org/eclipse/elk/) (for building pure Java projects that use ELK)



## Details

This is mainly a bugfix release. See GitHub for the full [list of resolved issues](https://github.com/eclipse-elk/elk/milestone/17?closed=1).


### New Features and Enhancements

* [#827](https://github.com/eclipse-elk/elk/pull/827): Added Greedy Model Order cycle breaker.
* [#842](https://github.com/eclipse-elk/elk/pull/842): Added a port label placement option to place the label always on the other same side.


### Bugfixes

* [#828](https://github.com/eclipse-elk/elk/pull/828): Node order violations are only counted for real nodes.
* [#829](https://github.com/eclipse-elk/elk/pull/829): Greedy Model Order cycle breaker is actually created if the option is set.
* [#825](https://github.com/eclipse-elk/elk/pull/825), [#826](https://github.com/eclipse-elk/elk/issues/826): Force model order option turns off when turned off.
* [#842](https://github.com/eclipse-elk/elk/pull/842): Vertical and horizontal port label spacing work for negative values.
