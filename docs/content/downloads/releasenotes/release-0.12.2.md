---
title: "0.12.2"
menu:
  main:
    parent: "ReleaseNotes"
    weight: -122
---

* [Release log](https://projects.eclipse.org/projects/modeling.elk/releases/0.12.2)
* [Documentation](https://download.eclipse.org/elk/updates/releases/0.12.2/elk-0.12.2-docs.zip)
* [Update site](https://download.eclipse.org/elk/updates/releases/0.12.2/)
* [Zipped update site](https://download.eclipse.org/elk/updates/releases/0.12.2/elk-0.12.2.zip) (for offline use)
* [Maven central](https://repo.maven.apache.org/maven2/org/eclipse/elk/) (for building pure Java projects that use ELK)



## Details

This release is a minor patch containing some bugfixes as well as some cleanup and upgrades under the hood. See GitHub for the full [list of resolved issues](https://github.com/eclipse-elk/elk/milestone/26?closed=1).

### Bugfixes
- [#1204](https://github.com/eclipse-elk/elk/issues/1204), [#1234](https://github.com/eclipse-elk/elk/pull/1234): Layer and position IDs should now be returned in the JSON result, when enabling generatePositionAndLayerIds.
- [#1211](https://github.com/eclipse-elk/elk/issues/1211), [#1242](https://github.com/eclipse-elk/elk/pull/1242): Shift coincident SPOrE node centres with a seeded random generator
- [#1214](https://github.com/eclipse-elk/elk/issues/1214), [#1215](https://github.com/eclipse-elk/elk/pull/1215): rectpacking: Handle more than two equal candidates in the area approximation
- [#1236](https://github.com/eclipse-elk/elk/issues/1236), [#1237](https://github.com/eclipse-elk/elk/pull/1237): Keep the interactive rectpacking tie order with a consistent comparator
- [#1216](https://github.com/eclipse-elk/elk/pull/1216), [#1217](https://github.com/eclipse-elk/elk/pull/1217): Fixes to EndLabelSorter
- [#1213](https://github.com/eclipse-elk/elk/pull/1213): Drop a repeated assignment in InLayerConstraintProcessorTest
- [#1210](https://github.com/eclipse-elk/elk/pull/1210): Read the port influence when choosing the model-order crossing counter
- [#1218](https://github.com/eclipse-elk/elk/issues/1218), [#1219](https://github.com/eclipse-elk/elk/pull/1219): Refuse BF/DF model order layerings with edges that do not point forward
- [#1220](https://github.com/eclipse-elk/elk/issues/1220), [#1221](https://github.com/eclipse-elk/elk/pull/1221): Ignore in-layer constraints whose target is not in the node's layer
- [#1222](https://github.com/eclipse-elk/elk/issues/1222), [#1223](https://github.com/eclipse-elk/elk/pull/1223), [#1079](https://github.com/eclipse-elk/elk/issues/1079), [#1245](https://github.com/eclipse-elk/elk/pull/1245): Fixes to graph wrapping. Graph wrapping can now handle when there are edges with labels on them
- [#1228](https://github.com/eclipse-elk/elk/issues/1228), [#1229](https://github.com/eclipse-elk/elk/pull/1229): Keep the edges of spline hyperedge segments in a deterministic order
- [#1224](https://github.com/eclipse-elk/elk/issues/1224), [#1225](https://github.com/eclipse-elk/elk/pull/1225): Read the rightmost port's side for four-sided self loop labels
- [#1232](https://github.com/eclipse-elk/elk/issues/1232), [#1233](https://github.com/eclipse-elk/elk/pull/1233): Fix inside self-loop label coordinates
- [#1226](https://github.com/eclipse-elk/elk/issues/1226), [#1227](https://github.com/eclipse-elk/elk/pull/1227): Let PortRestorer split self loop ports at the first non-negative flow
- [#1230](https://github.com/eclipse-elk/elk/issues/1230), [#1231](https://github.com/eclipse-elk/elk/pull/1231): Fix fixed port order for inside self-loops in compound layouts
- [#1246](https://github.com/eclipse-elk/elk/issues/1246), [#1247](https://github.com/eclipse-elk/elk/pull/1247): Report all overlaps of a node inserted into the overlap scanline
