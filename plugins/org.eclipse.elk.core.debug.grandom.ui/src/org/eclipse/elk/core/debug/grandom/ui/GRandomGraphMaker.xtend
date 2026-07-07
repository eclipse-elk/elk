/*******************************************************************************
 * Copyright (c) 2016, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.grandom.ui

import java.io.IOException
import java.util.ArrayList
import java.util.Collections
import java.util.List
import org.eclipse.core.resources.IFile
import org.eclipse.core.resources.IProject
import org.eclipse.core.runtime.CoreException
import org.eclipse.elk.core.debug.grandom.gRandom.Configuration
import org.eclipse.elk.core.debug.grandom.gRandom.Formats
import org.eclipse.elk.core.debug.grandom.gRandom.RandGraph
import org.eclipse.elk.core.debug.grandom.generators.RandomGraphGenerator
import org.eclipse.elk.graph.ElkNode
import org.eclipse.emf.common.util.URI
import org.eclipse.emf.ecore.resource.impl.ResourceSetImpl

class GRandomGraphMaker {

    val Iterable<Configuration> configs
    val DEFAULT_NAME = "random"
    val DEFAULT_FORMAT = Formats.ELKT.toString

    new(RandGraph rdg) {
        configs = rdg.configs
    }

    // here the generation of the graphs is triggered and the graphs are stored
    def gen(IProject project) {
        val generator = new RandomGraphGenerator(null)
        for (config : configs) {
            val filename = if(config.filename.exists) config.filename else DEFAULT_NAME
            val format = if(config.format.exists) config.format else DEFAULT_FORMAT
            val graphs = generator.generate(config)
            var fileNum = 0
            for (var i = 0; i < graphs.size; i++) {
                while (project.getFile(filename + fileNum + '.' + format).exists) {
                    fileNum++
                }
                val f = project.getFile(filename + fileNum + '.' + format)
                serialize(graphs.get(i), f)
            }
        }
    }

    //here the graphs are generated and returned instead of stored in a file
    def List<ElkNode> loadGraph() {
        val graphs = new ArrayList
        
        val generator = new RandomGraphGenerator(null)
        for (config : configs) {
            graphs.addAll(generator.generate(config))
        }
        return graphs
    }

    private def serialize(ElkNode graph, IFile file) throws IOException, CoreException {
        val resourceSet = new ResourceSetImpl();
        val resource = resourceSet.createResource(URI.createURI(file.getLocationURI().toString()));
        resource.getContents().add(graph);
        resource.save(Collections.EMPTY_MAP);
        file.refreshLocal(1, null);
    }

    private def exists(Object o) {
        return o !== null
    }
}
