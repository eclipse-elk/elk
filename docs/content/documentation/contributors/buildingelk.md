---
title: "Building ELK"
menu:
  main:
    identifier: "BuildingELK"
    parent: "ELKContributors"
    weight: 30
---

ELK is built using [Apache Maven](https://maven.apache.org/) in conjunction with [Tycho](https://eclipse.org/tycho/) to tell Maven how to build Eclipse projects. There are two parts that can be built: the Eclipse Layout Kernel itself, and the metadata compiler used by the main ELK build. The remainder of this page assumes that you have opened a shell in the `build/` directory inside your clone of the ELK repository.

## Build Requirements

Building ELK currently requires Java 8 and _will not_ work with more recent Java versions. We use Maven 3.6, but other versions will probably work just as well.


## Building ELK

Execute Maven using the following command line (note that the command line is split in order to improve readability):

```bash
mvn
    --define elkMeta.repositoryUrl=http://build.eclipse.org/modeling/elk/maven/meta/nightly
    --define elk.metadata.documentation.outputPath=/ELK_FOLDER/docs
    clean
    package
```

The meta repository URL above works if you haven't touched the meta data compiler and are trying to build the current master branch. If you have made changes to the meta data compiler, the URL should be the same that you used when compiling the compiler (see below). If you're trying to build a release, use the meta data compiler URL noted on the release's downloads page.


## Running Unit Tests

The unit tests require our [models repository](https://github.com/eclipse/elk-models) to be checked out as well. If it is, the following command line will run the tests:

```bash
mvn
    --define elkMeta.repositoryUrl=http://build.eclipse.org/modeling/elk/maven/meta/nightly
    --define elk.metadata.documentation.outputPath=/ELK_FOLDER/docs
    --define tests.paths.elk-repo=/ELK_FOLDER
    --define tests.paths.models-repo=/ELK_MODELS_FOLDER
    clean
    integration-test
```

The `TEST_RESULTS_FOLDER` is a folder where tests can save data about tests that have failed for later inspection.


## Building the Metadata Compiler

Execute Maven using the following command line to build the metadata compiler locally (note that the command line is split in order to improve readability):

```bash
mvn
    -P elk-meta
    --define elkMeta.repositoryUrl=file:///home/user/elkmetarepository
    --define elkMeta.publishUrl=file:///home/user/elkmetarepository
    clean
    deploy
```

## Building the Website and Documentation

The documentation website is our source of documentation. The main build produces reference documentation which is part of the website from our metadata files. The documentation website itself is then built with [Hugo](https://gohugo.io/), a static website generator. The website is setup in a way to be browsable both online and offline.

To actually build the documentation website, go to the repository's `docs/` folder and simply execute the following command:

```bash
hugo
```

To edit the website, start a Hugo server instance in the repository's `docs/` folder:

```bash
hugo server
```

You can then direct your browser to [http://localhost:1313/elk/](http://localhost:1313/elk/). As you edit the website's files and save them, the browser's content is dynamically updated.
