# Preparing a Release

## Release Review

The Eclipse release process is described in more detail in the [Eclipse Project Handbook](https://www.eclipse.org/projects/handbook/#release). The review process is tracked in an [issue such as this](https://gitlab.eclipse.org/eclipsefdn/emo-team/emo/-/issues/258).

1. Add a new release on ELK's [project page](https://projects.eclipse.org/projects/modeling.elk).
1. For major releases:
    1. Update NOTICE.md file with relevant information regarding used third-party content.
    2. Run `mvn org.eclipse.dash:license-tool-plugin:license-check -Dmaven.repo.local=./mvnrepo clean -Dtycho.target.eager=true -Dmaven.test.skip=true -Ddash.summary=../DEPENDENCIES` and commit the updated DEPENDENCIES file.
    3. Have [PMC](mailto:modeling-pmc@eclipse.org) approve the release documentation.
    4. Schedule a release review. Release reviews run for a minimum of one week and conclude on the first and third Wednesdays of each month.


## The Release Branch

1. Ensure that all bundles to be built and released are in fact built and are part of the [update site](https://download.eclipse.org/elk/updates/).
1. Create a release branch `releases/VERSION`.
1. Remove `-SNAPSHOT` and `.qualifier` from any version numbers. This is necessary for Maven to push the build to the proper Maven Central release staging area. Be careful not to change anything in source code; this should mainly affect `category.xml`, `pom.xml`s, `feature.xml`s and `MANIFEST.MF`s.
1. Open `build/org.eclipse.elk.repository/category.xml` and update its description like this:
   
    ```xml
    <description name="Eclipse Layout Kernel (Release VERSION_NUMBER)" url="https://download.eclipse.org/elk/updates/releases/VERSION_NUMBER">
      Update site for the Eclipse Layout Kernel, version VERSION_NUMBER.
    </description>
    ```
1. Remove the `Website` stage's call to the `publish-website.sh` script from the release build's `Jenkinsfile` and configure the build variables in the [Jenkins configuration](https://ci.eclipse.org/elk/job/IntegrationNightly/):

   
    Variable              | New value
    --------------------- | ---------------------------------------------------------------------
    `BRANCH`              | `releases/VERSION`
    `VERSION`             | Well... the version number...
    `ELK_TARGET_DIR`      | `/home/data/httpd/download.eclipse.org/elk/updates/releases/VERSION/`
    
    Pre-releases should point to a separate directory. Otherwise people that already use the pre-release have to clean their p2 directory to avoid problems.
    
1. \[deprecated\] Update the _[ReleaseNightly](https://ci.eclipse.org/elk/job/ReleaseNightly/)_ build with the same default values for the build variables. Also remember to update the repository branch the build will check out. Run it.
1. Update the version numbers on `master`. Tycho can help:

    ```
    mvn org.eclipse.tycho:tycho-versions-plugin:set-version -DnewVersion=NEW_VERSION_NUMBER.qualifier
    ```


## Releasing to Maven Central

This is a summary of the information on [this page](https://central.sonatype.org/pages/releasing-the-deployment.html).

1. Login to [OSSRH](https://oss.sonatype.org/).
1. In _Staging Repositories_, look for the drop that contains the new ELK bundles.
1. Examine their content.
1. Once satisfied, _Close_ the bundle. This will cause Nexus to examine the bundle and check whether everything is in order. If it isn't, _Drop_ the bundle, fix any problems and build a new release.
1. If closing the bundle was successful, _Release_ the bundle and be done with it.


## Releasing to Update Site

1. ~~Run the _PromoteUpdateSite_ build with proper parameter values and run it.~~
1. While in Jenkins, disable the nightly release build.


## The Website

1. Add a new website page for the release notes. It will appear on the website the next time the nightly build is run, so depending on the timing running that job manually might be an option.


## GitHub Management

1. Close the release milestone [on GitHub](https://github.com/eclipse-elk/elk/milestones).
1. Create a new release [on GitHub](https://github.com/eclipse-elk/elk/releases) with a release tag named `vVERSION` (for instance, `v0.4.0`).
