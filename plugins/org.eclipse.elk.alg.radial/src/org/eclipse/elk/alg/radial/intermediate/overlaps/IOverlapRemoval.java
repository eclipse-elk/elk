/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.radial.intermediate.overlaps;

import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * Interface for overlap removal strategies.
 *
 */
public interface IOverlapRemoval {

    /**
     * Remove the existing overlaps from the graph.
     * 
     * @param graph
     *          The graph where all nodes are positioned but overlaps may appear.
     * @param progressMonitor
     *          The progress monitor
     */
    void removeOverlaps(ElkNode graph, IElkProgressMonitor progressMonitor);
}
