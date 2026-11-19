---
title: "0.10.0"
menu:
  main:
    parent: "ReleaseNotes"
    weight: -100
---

* [Release log](https://projects.eclipse.org/projects/modeling.elk/releases/0.10.0)
* [Documentation](https://download.eclipse.org/elk/updates/releases/0.10.0/elk-0.10.0-docs.zip)
* [Update site](https://download.eclipse.org/elk/updates/releases/0.10.0/)
* [Zipped update site](https://download.eclipse.org/elk/updates/releases/0.10.0/elk-0.10.0.zip) (for offline use)
* [Maven central](https://repo.maven.apache.org/maven2/org/eclipse/elk/) (for building pure Java projects that use ELK)



## Details

This release concentrates on bug fixes and enhancements concerning many details of the layout algorithms as well as some cleanup and upgrades under the hood. ELK 0.10.0 now requires Java 17. See GitHub for the full [list of resolved issues](https://github.com/eclipse-elk/elk/milestone/20?closed=1).

### New Features and Enhancements
- [#1071](https://github.com/eclipse-elk/elk/pull/1071), [#1012](https://github.com/eclipse-elk/elk/pull/1012), [#901](https://github.com/eclipse-elk/elk/pull/901): Provide new properties for JSON coordinate options.
- [#1095](https://github.com/eclipse-elk/elk/pull/1095), [#1100](https://github.com/eclipse-elk/elk/pull/1100): Add OrderBySize Rectpacking preprocessor.
- [#1097](https://github.com/eclipse-elk/elk/pull/1097), [#1089](https://github.com/eclipse-elk/elk/pull/1089): Top-down layout can now be used without setting a size approximator.
- [#1094](https://github.com/eclipse-elk/elk/pull/1094), [#1052](https://github.com/eclipse-elk/elk/pull/1052): Add a layer splitting intermediate processor for Layered.
- [#1082](https://github.com/eclipse-elk/elk/pull/1082): Add fuzziness property for label soft wrapping.
- [#1043](https://github.com/eclipse-elk/elk/pull/1043), [#1042](https://github.com/eclipse-elk/elk/pull/1042): Add timeout property to libavoid.

### Changes
- [#1130](https://github.com/eclipse-elk/elk/pull/1130), [#1125](https://github.com/eclipse-elk/elk/pull/1125): Moved from Java 11 to Java 17 and upgraded to Tycho 4.

### Bugfixes
- [#1109](https://github.com/eclipse-elk/elk/pull/1109): Use ELK's classloader, not root classloader to discover services.
- [#1104](https://github.com/eclipse-elk/elk/pull/1104), [#1001](https://github.com/eclipse-elk/elk/pull/1001): Remove old unused edge sections in JSON graphs.
- [#1081](https://github.com/eclipse-elk/elk/pull/1081), [#953](https://github.com/eclipse-elk/elk/pull/953): Fix Self-loop label bug when using compaction.
- [#1060](https://github.com/eclipse-elk/elk/pull/1060), [#1059](https://github.com/eclipse-elk/elk/pull/1059): Fix an incorrect calculation of how to shift a mrtree layout.
- [#1067](https://github.com/eclipse-elk/elk/pull/1067): Fix always existing bug in model order port sorting.
- [#1054](https://github.com/eclipse-elk/elk/pull/1054): Make whitespace elimination in rectpacking disabled by default again.
- [#1046](https://github.com/eclipse-elk/elk/pull/1046): Include padding when setting node sizes in top-down layout.
- [#955](https://github.com/eclipse-elk/elk/pull/955): Consider node's size when computing available space in port placement.

### Cleanup
- [#1113](https://github.com/eclipse-elk/elk/pull/1113): Use target platform as file.
- [#1118](https://github.com/eclipse-elk/elk/pull/1118): Explicitly mark all Eclipse projects as UTF-8.
- [#1119](https://github.com/eclipse-elk/elk/pull/1119): Remove duplicate groupid and version.
- [#1111](https://github.com/eclipse-elk/elk/pull/1111): Keep whitespace in files to prevent unintended changes in PRs.
- [#1112](https://github.com/eclipse-elk/elk/pull/1112), [#1115](https://github.com/eclipse-elk/elk/pull/1115), [#1116](https://github.com/eclipse-elk/elk/pull/1116): Maven build cleanup.
- [#1117](https://github.com/eclipse-elk/elk/pull/1117): Add missing module names.
- [#1103](https://github.com/eclipse-elk/elk/pull/1103): Integrate license-check-plugin.
- [#1075](https://github.com/eclipse-elk/elk/pull/1075): Remove unused imports.
- [#1041](https://github.com/eclipse-elk/elk/pull/1041): Remove checkstyle natures from all projects.
