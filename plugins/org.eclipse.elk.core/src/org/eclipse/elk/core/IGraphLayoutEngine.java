/*******************************************************************************
 * Copyright (c) 2008, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core;

import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * A graph layout engine is able to perform automatic layout on a graph or parts of it.
 * 
 * @author swe
 */
public interface IGraphLayoutEngine {

    /**
     * Performs layout on the given layout graph.
     * 
     * @param layoutGraph
     *            the top-level node of the graph to be laid out
     * @param progressMonitor
     *            monitor to which progress of the layout algorithms is reported
     * @throws UnsupportedGraphException
     *             if the given graph is not supported by this algorithm
     * @throws UnsupportedConfigurationException
     *             if the layout configuration included in the graph is inconsistent or incompatible
     */
    void layout(ElkNode layoutGraph, IElkProgressMonitor progressMonitor);

}
