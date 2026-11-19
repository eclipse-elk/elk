/*******************************************************************************
 * Copyright (c) 2013, 2018 Kiel University and others.
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
package org.eclipse.elk.core.debug.grandom.generators;

import java.util.Random;

import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.MapPropertyHolder;
import org.eclipse.elk.graph.properties.Property;

/**
 * Property holder for random graph generator options.
 *
 * @author msp
 */
public class GeneratorOptions extends MapPropertyHolder {

    /** the serial version UID. */
    private static final long serialVersionUID = -917483559995737504L;

    /** the available file formats. */
    public static enum FileFormat {
        /** XMI-based format (kgx). */
        XMI,
        /** Xtext-based format (kgt). */
        XTEXT
    }

    /** the available graph types. */
    public static enum GraphType {
        /** custom graph. */
        CUSTOM,
        /** a tree. */
        TREE,
        /** a biconnected graph. */
        BICONNECTED,
        /** a triconnected graph. */
        TRICONNECTED,
        /** an acyclic graph without transitive edges. */
        ACYCLIC_NO_TRANSITIVE_EDGES,
        /** a bipartite graph. */
        BIPARTITE;
    }

    /** the possible ways to determine edges. */
    public enum EdgeDetermination {
        /** absolute number of edges in the graph. */
        ABSOLUTE,
        /** number of outgoing edges per node. */
        OUTGOING,
        /** relative number of edges (relative to n). */
        RELATIVE,
        /** density (relative to n^2). */
        DENSITY
    }

    //~~~~~~~~~~~~~~~~ Options for file creation

    /** the preference key for the number of graphs. */
    public static final IProperty<Integer> NUMBER_OF_GRAPHS = new Property<Integer>(
            "randomWizard.numberOfGraphs", 1);
    /** the preference key for the filename. */
    public static final IProperty<String> FILE_NAME = new Property<String>(
            "randomWizard.filename", "random.kgt");
    /** the preference key for the file format. */
    public static final IProperty<FileFormat> FILE_FORMAT = new Property<FileFormat>(
            "randomWizard.fileFormat", FileFormat.XTEXT);

    //~~~~~~~~~~~~~~~~ All other options depend on which GRAPH_TYPE has been selected

    /** option for the graph type. */
    public static final IProperty<GraphType> GRAPH_TYPE = new Property<GraphType>("basic.graphType",
            GraphType.CUSTOM);

    //~~~~~~~~~~~~~~~~ Options for all graph types

    /** option for the random number of nodes. */
    public static final IProperty<RandVal> NUMBER_OF_NODES = new Property<RandVal>(
            "basic.numNodes", RandVal.exact(10));
    /** option for specifying how to determine edges. */
    public static final IProperty<EdgeDetermination> EDGE_DETERMINATION =
            new Property<EdgeDetermination>("basic.edgeDetermination",
                    EdgeDetermination.ABSOLUTE);
    /** option for the absolute number of edges. */
    public static final IProperty<RandVal> EDGES_ABSOLUTE = new Property<RandVal>(
            "basic.numberOfEdges", RandVal.exact(20));
    /** option for the standard deviation in the absolute number of edges. */
    public static final IProperty<RandVal> RELATIVE_EDGES = new Property<RandVal>(
            "basic.relEdges", RandVal.allNil());
    /** option for standard deviation in the relative number of edges. */
    public static final IProperty<RandVal> DENSITY = new Property<RandVal>(
            "basic.density", RandVal.allNil());
    /** option that enables hierarchical graphs. */
    public static final IProperty<Boolean> ENABLE_HIERARCHY = new Property<Boolean>(
            "basic.enableHierarchy", false);
    /** option for the chance of creating a compound node. */
    public static final IProperty<Float> HIERARCHY_CHANCE = new Property<Float>(
            "basic.hierarchyChance", 0.05f, 0.0f, 1.0f);
    /** option for the maximum hierarchy level. */
    public static final IProperty<RandVal> MAX_HIERARCHY_LEVEL = new Property<RandVal>(
            "basic.maxHierarchyLevel", RandVal.exact(3));
    /** option for the factor to calculate the number of nodes in a compound node. */
    public static final IProperty<Float> HIERARCHY_NODES_FACTOR = new Property<Float>(
            "basic.hierarchyNodesFactor", 0.5f, 0.0f);
    /** option for the chance of creating a hypernode. */
    public static final IProperty<Float> HYPERNODE_CHANCE = new Property<Float>(
            "basic.hypernodeChance", 0.0f, 0.0f, 1.0f);
    /** option for using ports to connect nodes. */
    public static final IProperty<Boolean> ENABLE_PORTS = new Property<Boolean>(
            "basic.ports", false);
    /** option for the chance of edges to use already existing ports. */
    public static final IProperty<RandVal> USE_EXISTING_PORTS_CHANCE = new Property<RandVal>(
            "basic.useExistingPortsChance", RandVal.exact(0.3));
    /** option for allowing cross-hierarchy edges. */
    public static final IProperty<Boolean> CROSS_HIERARCHY_EDGES = new Property<Boolean>(
            "basic.crossHierarchyEdges", false);
    /** option for using a time-based randomization seed. */
    public static final IProperty<Boolean> TIME_BASED_RANDOMIZATION = new Property<Boolean>(
            "basic.timeBasedRandomization", true);
    /** option for the fixed randomization seed value. */
    public static final IProperty<Integer> RANDOMIZATION_SEED = new Property<Integer>(
            "basic.randomizationSeed", 0);
    /** option for setting the size of nodes. */
    public static final IProperty<Boolean> SET_NODE_SIZE = new Property<Boolean>(
            "basic.setNodeSize", true);
    /** option for the random width of nodes. */
    public static final IProperty<RandVal> NODE_WIDTH = new Property<RandVal>(
            "basic.minNodeWidth", RandVal.exact(30));
    /** option for the random height of nodes. */
    public static final IProperty<RandVal> NODE_HEIGHT = new Property<RandVal>(
            "basic.minNodeHeight", RandVal.exact(30));
    /** option for creating node labels. */
    public static final IProperty<Boolean> CREATE_NODE_LABELS = new Property<Boolean>(
            "basic.createNodeLabels", false);
    /** option for setting the size of ports. */
    public static final IProperty<Boolean> SET_PORT_SIZE = new Property<Boolean>(
            "basic.setPortSize", true);
    /** option for setting the size of ports. */
    public static final IProperty<RandVal> PORT_WIDTH = new Property<RandVal>(
            "basic.portWidth", RandVal.exact(4d));
    /** option for setting the size of ports. */
    public static final IProperty<RandVal> PORT_HEIGHT = new Property<RandVal>(
            "basic.portHeight", RandVal.exact(4d));
    /** option for setting the size of ports. */
    public static final IProperty<Boolean> CREATE_PORT_LABELS = new Property<Boolean>(
            "basic.createPortLabels", false);

    //~~~~~~~~~~~~~~~~ Layout options

    /** option for port constraints. */
    public static final IProperty<PortConstraints> PORT_CONSTRAINTS = new Property<PortConstraints>(
            "layout.portConstraints", PortConstraints.UNDEFINED);
    /** option for relative probability of incoming edges on the north side. */
    public static final IProperty<Integer> INCOMING_NORTH_SIDE = new Property<Integer>(
            "layout.incomingNorthSide", 10);
    /** option for relative probability of incoming edges on the east side. */
    public static final IProperty<Integer> INCOMING_EAST_SIDE = new Property<Integer>(
            "layout.incomingEastSide", 5);
    /** option for relative probability of incoming edges on the south side. */
    public static final IProperty<Integer> INCOMING_SOUTH_SIDE = new Property<Integer>(
            "layout.incomingSouthSide", 10);
    /** option for relative probability of incoming edges on the west side. */
    public static final IProperty<Integer> INCOMING_WEST_SIDE = new Property<Integer>(
            "layout.incomingWestSide", 75);
    /** option for relative probability of outgoing edges on the north side. */
    public static final IProperty<Integer> OUTGOING_NORTH_SIDE = new Property<Integer>(
            "layout.outgoingNorthSide", 10);
    /** option for relative probability of outgoing edges on the east side. */
    public static final IProperty<Integer> OUTGOING_EAST_SIDE = new Property<Integer>(
            "layout.outgoingEastSide", 75);
    /** option for relative probability of outgoing edges on the south side. */
    public static final IProperty<Integer> OUTGOING_SOUTH_SIDE = new Property<Integer>(
            "layout.outgoingSouthSide", 10);
    /** option for relative probability of outgoing edges on the west side. */
    public static final IProperty<Integer> OUTGOING_WEST_SIDE = new Property<Integer>(
            "layout.outgoingWestSide", 5);

    //~~~~~~~~~~~~~~~~ Options for GRAPH_TYPE CUSTOM

    /** option for the minimum number of outgoing edges. */
    public static final IProperty<RandVal> OUTGOING_EDGES = new Property<RandVal>(
            "basic.minOutgoingEdges", RandVal.allNil());
    /** option for allowing self-loops. */
    public static final IProperty<Boolean> SELF_LOOPS = new Property<Boolean>(
            "basic.selfLoops", false);
    /** option for allowing multi-edges. */
    public static final IProperty<Boolean> MULTI_EDGES = new Property<Boolean>(
            "basic.multiEdges", true);
    /** option for allowing cycles. */
    public static final IProperty<Boolean> CYCLES = new Property<Boolean>(
            "basic.cycles", true);
    /** option for allowing isolated nodes. */
    public static final IProperty<Boolean> ISOLATED_NODES = new Property<Boolean>(
            "basic.isolatedNodes", true);
    /** option for generating random edge labels. */
    public static final Property<Boolean> EDGE_LABELS = new Property<Boolean>(
            "basic.edgeLabels", false);

    //~~~~~~~~~~~~~~~~ Options for GRAPH_TYPE TREE

    /** option for the maximum tree width ('0' meaning unlimited). */
    public static final IProperty<Integer> MAX_WIDTH = new Property<Integer>("basic.maxWidth", 0);
    /** option for the maximum degree ('0' meaning unlimited). */
    public static final IProperty<Integer> MAX_DEGREE = new Property<Integer>("basic.maxDegree", 0);

    //~~~~~~~~~~~~~~~~ Options for GRAPH_TYPE ACYCLIC_NO_TRANSITIVE_EDGES

    /** option for planarity. */
    public static final IProperty<Boolean> PLANAR = new Property<Boolean>("basic.planar", false);

    //~~~~~~~~~~~~~~~~ Options for GRAPH_TYPE BIPARTITE

    /** option for minimal fraction of nodes in second partition set. */
    public static final IProperty<RandVal> PARTITION_FRAC = new Property<RandVal>(
            "basic.minPartitionFraction", RandVal.exact(0.4f));

    // ~~~~~~~~~~~~~~~~ Options for GRAPH_TYPE HIERARCHICAL
    public static final IProperty<Boolean> SMALL_HIERARCHY = new Property<Boolean>(
            "basic.SMALL_HIERARCHY", false);
    /** Number of hierarchical nodes per graph. */
    public static final IProperty<RandVal> NUMBER_HIERARCHICAL_NODES = new Property<RandVal>(
            "basic.HIERARCH_NODE", RandVal.allNil());
    //TODO-alan rename.
    public static final IProperty<RandVal> CROSS_HIER = new Property<RandVal>(
            "basic.CROSS_HIER", RandVal.allNil());
    public static final IProperty<RandVal> EXACT_RELATIVE_HIER = new Property<RandVal>(
            "basic.EXACT_RELATIVE_HIER", null);  

   

    public static class RandVal {
        private double min;
        private double max;
        private double mean;
        private double stddv;
        private double exact;
        private boolean useMinMax;
        private boolean useExact;

        private RandVal(final double min, final double max, final double mean, final double stddv, final double exact,
                final boolean useMinMax,
                final boolean useExact) {
            this.min = min;
            this.max = max;
            this.mean = mean;
            this.stddv = stddv;
            this.exact = exact;
            this.useMinMax = useMinMax;
            this.useExact = useExact;
        }

        public static RandVal allNil() {
            return new RandVal(0, 0, 0, 0, 0, false, false);
        }

        public static RandVal minMax(final double min, final double max) {
            return new RandVal(min, max, 0, 0, 0, true, false);
        }

        public static RandVal exact(final double exact) {
            return new RandVal(0, 0, 0, 0, exact, false, true);
        }

        public static RandVal gaussian(final double mean, final double stddv) {
            return new RandVal(0, 0, mean, stddv, 0, false, false);
        }

        public double val(final Random r) {
            double result = 0;
            if (useExact) {
                result = exact;
            } else if (useMinMax) {
                assert max >= min;
                result = r.nextDouble() * (max - min) + min;
            } else {
                result = r.nextGaussian() * stddv + mean;
            }
            return result < 0 ? 0 : result;
        }

        public int intVal(final Random r) {
            return (int) Math.round(val(r));
        }

        public int defaultInt() {
            return (int) Math.round(exact);
        }

        public float floatVal(final Random r) {
            return (float) val(r);
        }

        public void setMean(final double d) {
            mean = d;
        }

        public void setStddv(final double d) {
            stddv = d;
        }

        public void setMin(final int m) {
            min = m;
        }

        public void setMax(final int m) {
            max = m;
        }

        public double min() {
            return min;
        }

        public double max() {
            return max;
        }

        public double defaultVal() {
            return exact;
        }
    }

}
