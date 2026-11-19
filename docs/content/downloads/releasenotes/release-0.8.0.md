---
title: "0.8.0"
menu:
  main:
    parent: "ReleaseNotes"
    weight: -80
---

* [Release log](https://projects.eclipse.org/projects/modeling.elk/releases/0.8.0)
* [Documentation](https://download.eclipse.org/elk/updates/releases/0.8.0/elk-0.8.0-docs.zip)
* [Update site](https://download.eclipse.org/elk/updates/releases/0.8.0/)
* [Zipped update site](https://download.eclipse.org/elk/updates/releases/0.8.0/elk-0.8.0.zip) (for offline use)
* [Maven central](https://repo.maven.apache.org/maven2/org/eclipse/elk/) (for building pure Java projects that use ELK)



## Details

This is mainly a bugfix release. See GitHub for the full [list of resolved issues](https://github.com/eclipse-elk/elk/milestone/13).


### New Features and Enhancements

* [#672](https://github.com/eclipse-elk/elk/issues/672), [#674](https://github.com/eclipse-elk/elk/pull/674), [#675](https://github.com/eclipse-elk/elk/issues/675), [#677](https://github.com/eclipse-elk/elk/pull/677): Build systems is simplified.
* [#690](https://github.com/eclipse-elk/elk/issues/690), [#691](https://github.com/eclipse-elk/elk/issues/691): Improved documentation of position and layer choice constrains.
* [#695](https://github.com/eclipse-elk/elk/issues/695), [#698](https://github.com/eclipse-elk/elk/pull/698): Support node micro layout with further layout algorithms.
*  [#688](https://github.com/eclipse-elk/elk/issues/688), [#711](https://github.com/eclipse-elk/elk/pull/711): Better documentation for content alignment.
* [#722](https://github.com/eclipse-elk/elk/issues/722): Migrated to new build server.
* [#717](https://github.com/eclipse-elk/elk/pull/717): Model order: Property to weight model order node or port violations against edge crossings during crossing minimization. This also renames `considerModelOrder` to `considerModelOrder.strategy`.
* [#759](https://github.com/eclipse-elk/elk/pull/759): Model order: Added cycle breaker that enforces model order (but not against layerConstraints).
* [#815](https://github.com/eclipse-elk/elk/pull/815): Model order: Added option to enforce node order that existed before crossing minimization. This is to be used together with `considerModelOrder.strategy` :`NODES_AND_EDGES`
* [#816](https://github.com/eclipse-elk/elk/pull/816): Model order: Added property for nodes to signal that these should not get a model order and should be handled as dummy nodes.
* [#676](https://github.com/eclipse-elk/elk/issues/676), [#789](https://github.com/eclipse-elk/elk/pull/789): ServiceLoader can now be used with different class loaders.
* [#795](https://github.com/eclipse-elk/elk/pull/795): Added a crossing minimizer that does nothing instead of crossing minimization. This is to be used together to enforce the model order.
* [#804](https://github.com/eclipse-elk/elk/pull/804): Added fixed graph size support for the layered algorithm.
* [#780](https://github.com/eclipse-elk/elk/issues/780), [#802](https://github.com/eclipse-elk/elk/pull/802): Added option to generate position and layer ids.
* [#335](https://github.com/eclipse-elk/elk/issues/335), [#803](https://github.com/eclipse-elk/elk/pull/803): Added spacing documentation overview page to the website.
* [#819](https://github.com/eclipse-elk/elk/pull/819), [#822](https://github.com/eclipse-elk/elk/pull/822), [#823](https://github.com/eclipse-elk/elk/pull/823): Added option to order components with external ports not by their port connections but truly by model order.

### Changes

- [#717](https://github.com/eclipse-elk/elk/pull/717): Renamed `considerModelOrder` to `considerModelOrder.strategy`.
- [#697](https://github.com/eclipse-elk/elk/pull/697): Changed default node label stacking direction for UNDEFINED direction from horizontal to vertical.
- [#757](https://github.com/eclipse-elk/elk/issues/757), [#761](https://github.com/eclipse-elk/elk/pull/761), [#720](https://github.com/eclipse-elk/elk/issues/720), [#730](https://github.com/eclipse-elk/elk/pull/730): Bump guava version and remove upper bound.
- [#774](https://github.com/eclipse-elk/elk/pull/774): Update Eclipse download URL from http://build.eclipse.org/modeling/elk/updates to https://download.eclipse.org/elk/updates.
- [#651](https://github.com/eclipse-elk/elk/issues/651), [#790](https://github.com/eclipse-elk/elk/pull/790): Report invalid hierarchical crossing minimization instead of fixing it.
- [#791](https://github.com/eclipse-elk/elk/pull/791), [#714](https://github.com/eclipse-elk/elk/issues/714): Hierarchical edge orientation now uses the same in/out degree mechanism as the non-hierarchical layered algorithm.
- [#766](https://github.com/eclipse-elk/elk/issues/766), [#801](https://github.com/eclipse-elk/elk/pull/801), Hierarchical port dummies are considered to take no size, same as other port dummies.
- [#788](https://github.com/eclipse-elk/elk/pull/788): Individual spacings for vertical and horizontal label port spacing. `spacing.portLabel` is split into `spacing.labelPortHorizontal` and `spacing.labelPortVertical`.
- [#817](https://github.com/eclipse-elk/elk/pull/817), [#818](https://github.com/eclipse-elk/elk/pull/818): Separate connected components are no longer implicitly ordered my minimal model order of the component if a model order strategy is set but require a separate option.

### Removal

- [#526](https://github.com/eclipse-elk/elk/issues/526), [#760](https://github.com/eclipse-elk/elk/pull/760): Removed `layoutProvider` extension point.

### Bugfixes

- [#679](https://github.com/eclipse-elk/elk/pull/679): Prevent repeated registration of layout options during PlainJavaInitialization.

* [#706](https://github.com/eclipse-elk/elk/pull/706): Workaround for broken eclipse.urischeme dependency. Should be removed in next release once Eclipse dependencies are updated.
* [#718](https://github.com/eclipse-elk/elk/issues/718), [#721](https://github.com/eclipse-elk/elk/pull/721), [#727](https://github.com/eclipse-elk/elk/issues/727), [#728](https://github.com/eclipse-elk/elk/issues/728): Various fixes for the website.
* [#684](https://github.com/eclipse-elk/elk/issues/684), [#729](https://github.com/eclipse-elk/elk/pull/729): Adjusted melk file to show that a padding has to be applied to nodes.
* [#756](https://github.com/eclipse-elk/elk/pull/756): Model order is also set for hierarchical graphs.
* [#781](https://github.com/eclipse-elk/elk/issues/781), [#782](https://github.com/eclipse-elk/elk/pull/782): Clone mutable values instead of potentially sharing the same instance over multiple elements.
* [#784](https://github.com/eclipse-elk/elk/pull/784): Backward edges are correctly handled for `considerModelOrder`.
* [#775](https://github.com/eclipse-elk/elk/pull/775), [#776](https://github.com/eclipse-elk/elk/issues/776): Improved edge containment for JSON graphs.
* [#744](https://github.com/eclipse-elk/elk/issues/744), [#793](https://github.com/eclipse-elk/elk/pull/793): Radial layout now correctly calculates the radii im one radius is changed since nodes should not overlap.
* [#754](https://github.com/eclipse-elk/elk/issues/754), [#792](https://github.com/eclipse-elk/elk/pull/792): Dummy ports in hierarchical layered graphs get a default position of 0,0 to prevent an NPE.
* [#796](https://github.com/eclipse-elk/elk/pull/796): Graphviz no longer calls eclipse.ui if no platform is running.

### Cleanup

* [#694](https://github.com/eclipse-elk/elk/pull/694), [#704:](https://github.com/eclipse-elk/elk/pull/704) Cleanup.

### Known Bugs

- If `forceNodeModelOrder` is enabled once, it is still applied after it was disabled
