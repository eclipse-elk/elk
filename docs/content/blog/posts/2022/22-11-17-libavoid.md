---
title: "Edge Routing with Libavoid"
menu:
  main:
    identifier: "22-11-17-libavoid"
    parent: "2022"
    weight: 20
---
_By Miro Spönemann, November 17, 2022_

I'm happy to announce the availability of a new Maven module / Eclipse plug-in `org.eclipse.elk.alg.libavoid` that offers orthogonal edge routing with fixed nodes.
Credits for this work go to Ulf Rüegg, a former employee of the [Real-Time and Embedded Systems group](https://www.rtsys.informatik.uni-kiel.de/), where most of the development on ELK is happening.
My contribution is merely to revive this code and include it in the building and publishing process of the ELK project.

Routing edges without changing node positions has been a highly requested feature for several years.
This would enable users to freely move the nodes around and have the algorithm take care of the edges.
The existing layout algorithms such as _ELK Layered_ cannot do this because they have specific constraints for positioning the nodes.
Efforts to implement a standalone routing algorithm natively in Java have been started in the past, but they have not been brought to a working condition yet.
Therefore we decided to adopt the [libavoid](https://www.adaptagrams.org/documentation/libavoid.html) algorithm, which is part of the [Adaptagrams](http://www.adaptagrams.org) project.

There are two challenges for the integration of libavoid into ELK:

 * libavoid is written in C++ and built with [Autotools](https://www.gnu.org/software/automake/)
 * libavoid is licensed under the GNU Lesser General Public License (LGPL) version 2.1

## Integrating a C++ Library in ELK

Basically there are two ways to get a C++ based layout algorithm running with a Java application:

 * Use [SWIG](https://www.swig.org) two wrap the C++ library with JNI code to be accessed directly from the Java process
 * Compile the C++ library to an executable and communicate with it via standard I/O

After comparing these approaches, we found the one based on standard I/O to be more stable, so that's what we went for.
The [libavoid-server](https://github.com/TypeFox/libavoid-server) project reads graph descriptions and layout options from standard input, calls the libavoid router, and writes the resulting edge bend points to standard output.
The code is compiled with GitHub Actions for Linux, MacOS and Windows, and the generated executables are attached to [GitHub releases](https://github.com/TypeFox/libavoid-server/releases).
This project could be used independently of ELK to integrate libavoid with other kinds of applications.

## Integrating an LGPL Library in ELK

The Eclipse Foundation puts a strong emphasis for its projects to ensure that [intellectual property is handled properly](https://www.eclipse.org/projects/handbook/#ip).
This is important to avoid legal issues and to maintain a high level of trust from the users of an Eclipse project.

I went the formal path of requesting the usage of libavoid under LGPL license in an optional module of ELK.
The request was approved, as LGPL is generally considered to be compatible with the Eclipse Public License (EPL).

The new libavoid algorithm is offered as a separate feature named "Libavoid Integration for ELK" on the [ELK update site](/downloads.html), meaning that it is not automatically installed with the usual "Layout Algorithms" feature.
If you include this feature in your application, you need to comply with the [terms stated in the LGPL](https://github.com/mjwybrow/adaptagrams/blob/master/cola/LICENSE), like keeping the original source code accessible and distributing proper copyright and license notices.

## Using the Edge Router

The new layout algorithm has the identifier `org.eclipse.elk.alg.libavoid`. Its usage is basically the same as with other ELK layout algorithms, but there are a few specialties to consider:

 * Every node must have a position and a size (this algorithm won't move your nodes).
 * If your graph uses ports:
    * Every port must have a position on the border of its parent node.
    * `org.eclipse.elk.port.side` must be set to either `NORTH`, `EAST`, `SOUTH` or `WEST` for each port.

The libavoid algorithm also brings its own configuration options.
They are included in ELK's [reference documentation](/reference.html), but you can find documentation on the Doxygen pages of libavoid as well:

 * [RoutingOption](https://www.adaptagrams.org/documentation/namespaceAvoid.html#abc707ccbd6a0a7c29c124162c864ca05)
 * [RoutingParameter](https://www.adaptagrams.org/documentation/namespaceAvoid.html#a8a0154ae39129e7737d98e5a83daed19)

That's it for now. Please try this new integration and give us feedback on the [ELK repository](https://github.com/eclipse-elk/elk)!
