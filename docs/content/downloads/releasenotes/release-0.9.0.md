---
title: "0.9.0"
menu:
  main:
    parent: "ReleaseNotes"
    weight: -90
---

* [Release log](https://projects.eclipse.org/projects/modeling.elk/releases/0.9.0)
* [Documentation](https://download.eclipse.org/elk/updates/releases/0.9.0/elk-0.9.0-docs.zip)
* [Update site](https://download.eclipse.org/elk/updates/releases/0.9.0/)
* [Zipped update site](https://download.eclipse.org/elk/updates/releases/0.9.0/elk-0.9.0.zip) (for offline use)
* [Maven central](https://repo.maven.apache.org/maven2/org/eclipse/elk/) (for building pure Java projects that use ELK)



## Details

This is mainly a bugfix release. See GitHub for the full [list of resolved issues](https://github.com/eclipse-elk/elk/milestone/15?closed=1).


### New Features and Enhancements

* [#962](https://github.com/eclipse-elk/elk/pull/962), [#914](https://github.com/eclipse-elk/elk/pull/914): Added an experimental Depth-First and Breadth-First model order layerer.
* [#867](https://github.com/eclipse-elk/elk/pull/867), [#902](https://github.com/eclipse-elk/elk/pull/902/): Added a model order layering by node promotion.
* [#956](https://github.com/eclipse-elk/elk/pull/956), [#942](https://github.com/eclipse-elk/elk/pull/942), [#927](https://github.com/eclipse-elk/elk/pull/927), [#926](https://github.com/eclipse-elk/elk/pull/926), [#921](https://github.com/eclipse-elk/elk/pull/921), [#893](https://github.com/eclipse-elk/elk/pull/893), [#892](https://github.com/eclipse-elk/elk/pull/892), [#886](https://github.com/eclipse-elk/elk/pull/886): Added libavoid for standalone edge routing.
* [#945](https://github.com/eclipse-elk/elk/pull/945), [#941](https://github.com/eclipse-elk/elk/pull/941): Added option for initial rotation of radial layout.
* [#947](https://github.com/eclipse-elk/elk/issues/947), [#951](https://github.com/eclipse-elk/elk/pull/951), [#940](https://github.com/eclipse-elk/elk/pull/940), [#936](https://github.com/eclipse-elk/elk/pull/936), [\#932,](https://github.com/eclipse-elk/elk/pull/932) [#843](https://github.com/eclipse-elk/elk/pull/843): Added topdown layout and top-down rectangle packing algorithm.
* [#939](https://github.com/eclipse-elk/elk/pull/939), [#930](https://github.com/eclipse-elk/elk/issues/930): Added port model order.
* [#922](https://github.com/eclipse-elk/elk/pull/922): KVectors now have a rotation utility method.
* [#406](https://github.com/eclipse-elk/elk/issues/406), [#884](https://github.com/eclipse-elk/elk/pull/884): MrTree edge routing options + interactivity and model order.
* [#889](https://github.com/eclipse-elk/elk/pull/889): Added experimental bendpoints to force layout.
* [#877:](https://github.com/eclipse-elk/elk/pull/877) Added new model order strategy. 
* [#834](https://github.com/eclipse-elk/elk/pull/834), [#833](https://github.com/eclipse-elk/elk/issues/833), [#855](https://github.com/eclipse-elk/elk/pull/855), [#843](https://github.com/eclipse-elk/elk/pull/843):  Support for fixed graphs for mrtree, force, stress, radial, and rectpacking.

### Changes

- [#973](https://github.com/eclipse-elk/elk/pull/973), [#702](https://github.com/eclipse-elk/elk/issues/702), [#794](https://github.com/eclipse-elk/elk/issues/794), [#900](https://github.com/eclipse-elk/elk/pull/900), [#807](https://github.com/eclipse-elk/elk/issues/807): Moved from Java 8 to Java 11 while testing compliance with Java 17.
-  [#876](https://github.com/eclipse-elk/elk/pull/876): Renamed and reworked FORCE_MODEL_ORDER component ordering strategy to MODEL_ORDER and GROUP_MODEL_ORDER component ordering strategy.
- [#651](https://github.com/eclipse-elk/elk/issues/651), [#885](https://github.com/eclipse-elk/elk/pull/885): Report invalid hierarchical cross min configuration instead of fixing it.
- [#865](https://github.com/eclipse-elk/elk/pull/865), [#835](https://github.com/eclipse-elk/elk/pull/835): Restructure and improve rectpacking. This includes a lot of renaming. `optimizationGoal` is renamed to `widthApproximation.optimizationGoal`, `lastPlaceShift` is renamed to `widthApproximation.lastPlaceShift`, `onlyFirstIteration` is deleted, `rowCompaction` is now the default `packing.strategy` `COMPACTION`, `packing.compaction.rowHeightReevaluation` and `packing.compaction.iterations` are added to further improve the packing compaction, `expandNodes` and `expandToAspectRatio` have been reworked into `whiteSpaceElimination.strategy`s `EQUAL_BETWEEN_STRUCTURES` and `TO_ASPECT_RATIO`, `targetWidth` is now set via `widthApproximation.targetWidth`and is enabled by setting the `widthApproximation.strategy`, the `inNewRow` option can constrain the packing, and the `tryBox` option may be set to try whether the rectangles are stackable and in the case they are not use the simpler `box` layout algorithm.


### Bugfixes

* [#957](https://github.com/eclipse-elk/elk/pull/957): Corrected node node spacing usage for radial layout. The node node spacing is now halved and used as a padding around the nodes.
* [#946](https://github.com/eclipse-elk/elk/pull/946/), [#944](https://github.com/eclipse-elk/elk/issues/944): Fix hyper edge segment splitter to actually save the best area after finding it.
* [#850](https://github.com/eclipse-elk/elk/issues/850), [#938](https://github.com/eclipse-elk/elk/pull/938): Port labels no longer overlap with self loops.
* [#918](https://github.com/eclipse-elk/elk/issues/918), [#933](https://github.com/eclipse-elk/elk/pull/933): Fix port model order comparator.
* [#923](https://github.com/eclipse-elk/elk/pull/923): Wrapping strategies do now correctly work together with model order.
* [#906](https://github.com/eclipse-elk/elk/pull/906), [#905](https://github.com/eclipse-elk/elk/issues/905), [#913](https://github.com/eclipse-elk/elk/pull/913): Fix tail edge labels.
* [#869](https://github.com/eclipse-elk/elk/issues/869), [#911](https://github.com/eclipse-elk/elk/pull/911), [#868](https://github.com/eclipse-elk/elk/issues/868):  Model order bug fixes.
* [#870](https://github.com/eclipse-elk/elk/issues/870): Fixed NoSuchElementException in network simplex node placer.
* [#890](https://github.com/eclipse-elk/elk/pull/890): Correctly check whether the layered algorithm is set.
* [#852](https://github.com/eclipse-elk/elk/pull/852), [#841](https://github.com/eclipse-elk/elk/issues/841): Inline edge labels are now correctly centered (hopefully).
* [#887](https://github.com/eclipse-elk/elk/issues/887): Miscellaneous bug fixes.

### Cleanup

- [#924](https://github.com/eclipse-elk/elk/pull/924): Corrected pom license for maven central.



