/*******************************************************************************
 * Copyright (c) 2015, 2019 Kiel University and others.
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
package org.eclipse.elk.alg.layered;

import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.p4nodes.LinearSegmentsNodePlacer.LinearSegment;
import org.eclipse.elk.alg.layered.p5edges.orthogonal.HyperEdgeSegment;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.LoggedGraph;

/**
 * A utility class for debug output of ELK Layered. Currently it delegates to both {@link JsonDebugUtil} and
 * {@link DotDebugUtil}.
 */
public final class DebugUtil {

    private DebugUtil() {
    };
    
    /**
     * Logs a representation of the given graph in dot and JSON format.
     * 
     * @param monitor
     *            the progress monitor to log to.
     * @param lgraph
     *            the layered graph
     * @param slotIndex
     *            the slot before whose execution the graph is written.
     * @param name
     *            the name the slot before whose execution the graph is written.
     * @see {@link DotDebugUtil#createDebugGraph(LGraph, int)}
     * @see {@link JsonDebugUtil#createDebugGraph(LGraph, int)}
     */
    public static void logDebugGraph(final IElkProgressMonitor monitor, final LGraph lgraph, final int slotIndex,
            final String name) {
        
        String tag = slotIndex + "-" + name;
        monitor.logGraph(
                DotDebugUtil.createDebugGraph(lgraph),
                tag,
                LoggedGraph.Type.DOT);
        monitor.logGraph(
                JsonDebugUtil.createDebugGraph(lgraph),
                tag,
                LoggedGraph.Type.JSON);
    }

    /**
     * Logs a debug graph for the given linear segments and their dependencies.
     * 
     * @param monitor
     *            the progress monitor to log to.
     * @param layeredGraph
     *            the layered graph.
     * @param segmentList
     *            the list of linear segments.
     * @param outgoingList
     *            the list of successors for each linear segment.
     * @see {@link DotDebugUtil#createDebugGraph(LGraph, List, List)}
     * @see {@link JsonDebugUtil#createDebugGraph(LGraph, List, List)}
     */
    public static void logDebugGraph(final IElkProgressMonitor monitor, final LGraph layeredGraph,
            final List<LinearSegment> segmentList, final List<List<LinearSegment>> outgoingList) {

        String tag = "segment graph";
        monitor.logGraph(
                DotDebugUtil.createDebugGraph(layeredGraph, segmentList, outgoingList),
                tag,
                LoggedGraph.Type.DOT);
        monitor.logGraph(
                JsonDebugUtil.createDebugGraph(layeredGraph, segmentList, outgoingList),
                tag,
                LoggedGraph.Type.JSON);
    }

    /**
     * Logs a debug graph for the given list of hypernodes.
     * 
     * @param monitor
     *            the progress monitor to log to.
     * @param layeredGraph
     *            the layered graph
     * @param layerIndex
     *            the currently processed layer's index
     * @param hypernodes
     *            a list of hypernodes
     * @param debugPrefix
     *            prefix of debug output files
     * @param label
     *            a label to append to the output files
     * @see {@link DotDebugUtil#createDebugGraph(LGraph, int, List, String, String)}
     * @see {@link JsonDebugUtil#createDebugGraph(LGraph, int, List, String, String)}
     */
    public static void logDebugGraph(final IElkProgressMonitor monitor, final LGraph layeredGraph,
            final int layerIndex, final List<HyperEdgeSegment> hypernodes, final String debugPrefix,
            final String label) {

        String tag = debugPrefix + " - " + label;
        monitor.logGraph(
                DotDebugUtil.createDebugGraph(layeredGraph, hypernodes),
                tag,
                LoggedGraph.Type.DOT);
        monitor.logGraph(
                JsonDebugUtil.createDebugGraph(layeredGraph, hypernodes),
                tag,
                LoggedGraph.Type.JSON);
    }

}
