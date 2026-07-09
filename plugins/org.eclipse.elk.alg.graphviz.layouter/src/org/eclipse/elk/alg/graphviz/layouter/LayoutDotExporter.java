/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.alg.graphviz.layouter;

import java.util.List;

import org.eclipse.elk.alg.graphviz.dot.dot.Attribute;
import org.eclipse.elk.alg.graphviz.dot.dot.AttributeStatement;
import org.eclipse.elk.alg.graphviz.dot.dot.AttributeType;
import org.eclipse.elk.alg.graphviz.dot.dot.DotFactory;
import org.eclipse.elk.alg.graphviz.dot.dot.GraphvizModel;
import org.eclipse.elk.alg.graphviz.dot.dot.Statement;
import org.eclipse.elk.alg.graphviz.dot.transform.Attributes;
import org.eclipse.elk.alg.graphviz.dot.transform.Command;
import org.eclipse.elk.alg.graphviz.dot.transform.DotExporter;
import org.eclipse.elk.alg.graphviz.dot.transform.IDotTransformationData;
import org.eclipse.elk.alg.graphviz.dot.transform.NeatoModel;
import org.eclipse.elk.alg.graphviz.dot.transform.OverlapMode;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.EdgeLabelPlacement;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.options.HierarchyHandling;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;

import com.google.common.collect.Iterables;

/**
 * Specialized exporter for the Dot format that understands layout options.
 */
public class LayoutDotExporter extends DotExporter {
    
    @Override
    public void transform(final IDotTransformationData<ElkNode, GraphvizModel> transData) {
        ElkNode sourceGraph = transData.getSourceGraph();
        Command command = transData.getProperty(COMMAND);
        transData.setProperty(USE_EDGE_IDS, true);
        transData.setProperty(HIERARCHY,
                sourceGraph.getProperty(CoreOptions.HIERARCHY_HANDLING) == HierarchyHandling.INCLUDE_CHILDREN
                        && command == Command.DOT);
        transData.setProperty(TRANSFORM_NODE_LAYOUT, sourceGraph.getProperty(CoreOptions.INTERACTIVE));
        transData.setProperty(TRANSFORM_EDGE_LAYOUT, false);
        transData.setProperty(TRANSFORM_NODE_LABELS, false);
        transData.setProperty(ADAPT_PORT_POSITIONS, sourceGraph.getProperty(
                GraphvizMetaDataProvider.ADAPT_PORT_POSITIONS));
        super.transform(transData);
    }
    
    /** Base factor for setting {@link Attributes#SIMPLEX_LIMIT} from the iterations limit. */
    private static final float NSLIMIT_BASE = 100.0f;
    
    @Override
    protected void setGraphAttributes(final List<Statement> statements, final ElkNode parentNode,
            final IDotTransformationData<ElkNode, GraphvizModel> transData) {
        
        Command command = transData.getProperty(COMMAND);
        AttributeStatement graphAttrStatement = DotFactory.eINSTANCE.createAttributeStatement();
        graphAttrStatement.setType(AttributeType.GRAPH);
        List<Attribute> graphAttrs = graphAttrStatement.getAttributes();
        statements.add(graphAttrStatement);

        setGeneralNodeAttributes(statements);
        List<Attribute> edgeAttrs = setGeneralEdgeAttributes(statements);

        // set minimal spacing
        Double spacing = parentNode.getProperty(CoreOptions.SPACING_NODE_NODE);
        if (spacing == null || spacing < 0) {
            switch (command) {
            case CIRCO:
            case FDP:
            case NEATO:
                spacing = DEF_SPACING_LARGE;
                break;
            case TWOPI:
                spacing = DEF_SPACING_XLARGE;
                break;
            default:
                spacing = DEF_SPACING_SMALL;
            }
        }

        switch (command) {
        case DOT:
            graphAttrs.add(createAttribute(Attributes.NODESEP, spacing / DPI));
            double rankSepFactor = parentNode.getProperty(GraphvizMetaDataProvider.LAYER_SPACING_FACTOR);
            graphAttrs.add(createAttribute(Attributes.RANKSEP, rankSepFactor * spacing / DPI));
            // set layout direction
            switch (parentNode.getProperty(CoreOptions.DIRECTION)) {
            case UP:
                graphAttrs.add(createAttribute(Attributes.RANKDIR, "BT"));
                break;
            case LEFT:
                graphAttrs.add(createAttribute(Attributes.RANKDIR, "RL"));
                break;
            case RIGHT:
                graphAttrs.add(createAttribute(Attributes.RANKDIR, "LR"));
                break;
            default: // DOWN
                graphAttrs.add(createAttribute(Attributes.RANKDIR, "TB"));
            }
            // set iterations limit
            Double iterationsFactor = parentNode.getProperty(GraphvizMetaDataProvider.ITERATIONS_FACTOR);
            if (iterationsFactor != null && iterationsFactor > 0) {
                graphAttrs.add(createAttribute(Attributes.CROSSMIN_LIMIT, iterationsFactor));
                if (iterationsFactor < 1) {
                    double simplexLimit = iterationsFactor * NSLIMIT_BASE;
                    graphAttrs.add(createAttribute(Attributes.SIMPLEX_LIMIT, simplexLimit));
                }
            }
            // enable compound mode
            if (parentNode.getProperty(CoreOptions.HIERARCHY_HANDLING) == HierarchyHandling.INCLUDE_CHILDREN) {
                graphAttrs.add(createAttribute(Attributes.COMPOUND, "true"));
            }
            break;

        case TWOPI:
            graphAttrs.add(createAttribute(Attributes.RANKSEP, spacing / DPI));
            break;

        case CIRCO:
            graphAttrs.add(createAttribute(Attributes.MINDIST, spacing / DPI));
            break;

        case NEATO:
            edgeAttrs.add(createAttribute(Attributes.EDGELEN, spacing / DPI));
            // configure initial placement of nodes
            Integer seed = parentNode.getProperty(CoreOptions.RANDOM_SEED);
            if (seed == null) {
                seed = 1;
            } else if (seed == 0) {
                seed = -1;
            } else if (seed < 0) {
                seed = -seed;
            }
            graphAttrs.add(createAttribute(Attributes.START, "random" + seed));
            // set epsilon value
            Double epsilon = parentNode.getProperty(GraphvizMetaDataProvider.EPSILON);
            if (epsilon != null && epsilon > 0) {
                graphAttrs.add(createAttribute(Attributes.EPSILON, epsilon));
            }
            // set distance model
            NeatoModel model = parentNode.getProperty(GraphvizMetaDataProvider.NEATO_MODEL);
            if (model != NeatoModel.SHORTPATH) {
                graphAttrs.add(createAttribute(Attributes.NEATO_MODEL, model.literal()));
            }
            break;

        case FDP:
            graphAttrs.add(createAttribute(Attributes.SPRING_CONSTANT, spacing / DPI));
            break;
        }

        if (command == Command.NEATO || command == Command.FDP) {
            // set maximum number of iterations
            Integer maxiter = parentNode.getProperty(GraphvizMetaDataProvider.MAXITER);
            if (maxiter != null && maxiter > 0) {
                graphAttrs.add(createAttribute(Attributes.MAXITER, maxiter));
            }
        }

        if (command != Command.DOT) {
            // enable or disable node overlap avoidance
            OverlapMode mode = parentNode.getProperty(GraphvizMetaDataProvider.OVERLAP_MODE);
            if (mode != OverlapMode.NONE) {
                graphAttrs.add(createAttribute(Attributes.OVERLAP, mode.literal()));
                graphAttrs.add(createAttribute(Attributes.SEP, "\"+" + Math.round(spacing / 2)
                        + "\""));
            }
            // enable or disable connected component packing
            Boolean pack = parentNode.getProperty(CoreOptions.SEPARATE_CONNECTED_COMPONENTS);
            if (command == Command.TWOPI || pack != null && pack.booleanValue()) {
                graphAttrs.add(createAttribute(Attributes.PACK, spacing.intValue()));
            }
        }

        // configure edge routing
        EdgeRouting edgeRouting = parentNode.getProperty(CoreOptions.EDGE_ROUTING);
        String splineMode;
        switch (edgeRouting) {
        case POLYLINE:
            splineMode = "polyline";
            break;
        case ORTHOGONAL:
            splineMode = "ortho";
            break;
        default:
            splineMode = "spline";
            transData.setProperty(USE_SPLINES, true);
        }
        graphAttrs.add(createAttribute(Attributes.SPLINES, splineMode));

        // enable edge concentration
        if (parentNode.getProperty(GraphvizMetaDataProvider.CONCENTRATE)) {
            graphAttrs.add(createAttribute(Attributes.CONCENTRATE, "true"));
        }
    }
    
    @Override
    protected void setEdgeLabels(final ElkEdge kedge, final List<Attribute> attributes, final boolean isVertical) {
        super.setEdgeLabels(kedge, attributes, isVertical);
        // set label distance and angle
        if (Iterables.any(kedge.getLabels(), (ElkLabel label) -> {
            EdgeLabelPlacement elp = label.getProperty(CoreOptions.EDGE_LABELS_PLACEMENT);
            return elp == EdgeLabelPlacement.HEAD || elp == EdgeLabelPlacement.TAIL;
        })) {
            double distance = kedge.getProperty(GraphvizMetaDataProvider.LABEL_DISTANCE);
            if (distance >= 0.0) {
                attributes.add(createAttribute(Attributes.LABELDISTANCE, distance));
            }
            double angle = kedge.getProperty(GraphvizMetaDataProvider.LABEL_ANGLE);
            attributes.add(createAttribute(Attributes.LABELANGLE, angle));
        }
    }

}
