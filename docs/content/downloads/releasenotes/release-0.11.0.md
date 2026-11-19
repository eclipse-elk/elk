---
title: "0.11.0"
menu:
  main:
    parent: "ReleaseNotes"
    weight: -110
---

* [Release log](https://projects.eclipse.org/projects/modeling.elk/releases/0.11.0)
* [Documentation](https://download.eclipse.org/elk/updates/releases/0.11.0/elk-0.11.0-docs.zip)
* [Update site](https://download.eclipse.org/elk/updates/releases/0.11.0/)
* [Zipped update site](https://download.eclipse.org/elk/updates/releases/0.11.0/elk-0.11.0.zip) (for offline use)
* [Maven central](https://repo.maven.apache.org/maven2/org/eclipse/elk/) (for building pure Java projects that use ELK)



## Details

This release concentrates on bug fixes and enhancements concerning many details of the layout algorithms as well as some cleanup and upgrades under the hood. See GitHub for the full [list of resolved issues](https://github.com/eclipse-elk/elk/milestone/18?closed=1).

### New Features and Enhancements
- [#861](https://github.com/eclipse-elk/elk/issues/861), [#960](https://github.com/eclipse-elk/elk/pull/960): Add new algorithm 'VertiFlex' to layout trees with vertical position constraints for nodes.
- [#1147](https://github.com/eclipse-elk/elk/pull/1147), [#952](https://github.com/eclipse-elk/elk/issues/952): Add group model order strategies.
- [#1000](https://github.com/eclipse-elk/elk/issues/1000), [#1161](https://github.com/eclipse-elk/elk/pull/1161): CoffmanGrahamLayerer: Use node id as secondary order criterion.
- [#1145](https://github.com/eclipse-elk/elk/pull/1145), [#1163](https://github.com/eclipse-elk/elk/pull/1163): Add helper methods to log an LGraph.
- [#1096](https://github.com/eclipse-elk/elk/issues/1096), [#1129](https://github.com/eclipse-elk/elk/pull/1129): libavoidCancelRouting: Add a method that enables canceling of a started layout operation with libavoid.
- [#1144](https://github.com/eclipse-elk/elk/pull/1144): Add new size approximators for top-down layout.
- [#1107](https://github.com/eclipse-elk/elk/pull/1107): Add median heuristic for crossing minimization phase of layered.
- [#1150](https://github.com/eclipse-elk/elk/pull/1150): Add edge length minimizer heuristic for layer splitting.

### Changes
- [#866](https://github.com/eclipse-elk/elk/issues/866), [#1087](https://github.com/eclipse-elk/elk/pull/1087): Improve Box algorithm by considering the padding.

### Bugfixes
- [#659](https://github.com/eclipse-elk/elk/issues/659), [#1158](https://github.com/eclipse-elk/elk/pull/1158), [#1023](https://github.com/eclipse-elk/elk/issues/1023): Inline edge labels on self-loops with opposing sides.
- [#969](https://github.com/eclipse-elk/elk/issues/969), [#993](https://github.com/eclipse-elk/elk/issues/993), [#1160](https://github.com/eclipse-elk/elk/pull/1160): Layered: Make partitioning consider not directly connected partitions.
- [#1004](https://github.com/eclipse-elk/elk/issues/1004): Model order barycenter heuristic no longer breaks port order in hierachical graphs.
- [#968](https://github.com/eclipse-elk/elk/pull/968), [#1159](https://github.com/eclipse-elk/elk/pull/1159): Only reorder ports on parents if they changed.
- [#1140](https://github.com/eclipse-elk/elk/issues/1140), [#1142](https://github.com/eclipse-elk/elk/issues/1142): Recompute graph bounds after packing components in MrTree.
- [#1153](https://github.com/eclipse-elk/elk/issues/1153), [#1154](https://github.com/eclipse-elk/elk/issues/1154): Consider the parent label when having child nodes for size in Rectpacking.
- [#1152](https://github.com/eclipse-elk/elk/pull/1152): Fix incorrect node size calculation in topdown fallback case.
- [#1077](https://github.com/eclipse-elk/elk/issues/1077), [#1155](https://github.com/eclipse-elk/elk/issues/1155): rectpacking: Fix error in rectangle absorption.

### Website
- [#1133](https://github.com/eclipse-elk/elk/issues/1133), [#] : Add elkt format documentation to website.