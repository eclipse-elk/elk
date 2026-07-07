/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.radial.p1position;

import java.util.List;

import org.eclipse.elk.alg.radial.InternalProperties;
import org.eclipse.elk.alg.radial.RadialLayoutPhases;
import org.eclipse.elk.alg.radial.RadialUtil;
import org.eclipse.elk.alg.radial.intermediate.optimization.IEvaluation;
import org.eclipse.elk.alg.radial.options.RadialOptions;
import org.eclipse.elk.alg.radial.p1position.wedge.IAnnulusWedgeCriteria;
import org.eclipse.elk.alg.radial.sorting.IRadialSorter;
import org.eclipse.elk.core.alg.ILayoutPhase;
import org.eclipse.elk.core.alg.LayoutProcessorConfiguration;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * The initial radial layout is based on the algorithm of Peter Eades published in "Drawing free trees.", published by
 * International Institute for Advanced Study of Social Information Science, Fujitsu Limited in 1991. The radial
 * layouter takes a tree and places the nodes in radial order around the root. The nodes of the same tree level are
 * placed on the same radius. *
 */
public class EadesRadial implements ILayoutPhase<RadialLayoutPhases, ElkNode> {

    private static final int CIRCLE_DEGREES = 360;
    private static final double DEGREE_TO_RAD = Math.PI / 180;
    private double radius;
    private IRadialSorter sorter;
    private IAnnulusWedgeCriteria annulusWedgeCriteria;
    private IEvaluation optimizer;
    private ElkNode root;

    @Override
    public void process(final ElkNode graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Eades radial", 1);
        progressMonitor.logGraph(graph, "After");
        root = graph.getProperty(InternalProperties.ROOT_NODE);
        radius = graph.getProperty(RadialOptions.RADIUS);
        sorter = graph.getProperty(RadialOptions.SORTER).create();
        annulusWedgeCriteria = graph.getProperty(RadialOptions.WEDGE_CRITERIA).create();
        optimizer = graph.getProperty(RadialOptions.OPTIMIZATION_CRITERIA).create();
        translate(graph);
        progressMonitor.logGraph(graph, "After");
    }

    /**
     * Search for the best layout translation by looking at each degree.
     * 
     * @param root
     *            The root node of the graph
     */
    private void translate(final ElkNode graph) {
        double optimalOffset = 0;
        double optimalValue = Double.MAX_VALUE;

        if (optimizer != null) {
            for (int i = 0; i < CIRCLE_DEGREES; i++) {
                double offset = i * DEGREE_TO_RAD;
                positionNodes(root, 0, 0, 2 * Math.PI, offset);
                double translatedValue = optimizer.evaluate(root);
                // Take the first occurence of the minimum
                if (translatedValue < optimalValue) {
                    optimalOffset = offset;
                    optimalValue = translatedValue;
                }
            }
        }
        positionNodes(root, 0, 0, 2 * Math.PI, optimalOffset);
    }

    /**
     * Place a node in the center of a wedge and calculate the wedge for the next child.
     * 
     * @param node
     * @param currentRadius
     * @param minAlpha
     * @param maxAlpha
     * @param optimalOffset
     */
    private void positionNodes(final ElkNode node, final double currentRadius, final double minAlpha,
            final double maxAlpha, double optimalOffset) {
        double radOffest = optimalOffset;

        double alphaPoint = (minAlpha + maxAlpha) / 2 + radOffest;

        // x=r*sinθ, y=r*cosθ
        double xPos = currentRadius * Math.cos(alphaPoint);
        double yPos = currentRadius * Math.sin(alphaPoint);

        // shift the nodes, such that the center of each node is on the circle
        RadialUtil.centerNodesOnRadi(node, xPos, yPos);

        double numberOfLeafs = annulusWedgeCriteria.calculateWedgeSpace(node);
        double tau = 2 * Math.acos(currentRadius / currentRadius + radius);
        double s;
        double alpha;
        if (tau < maxAlpha - minAlpha) {
            s = tau / numberOfLeafs;
            alpha = (minAlpha + maxAlpha - tau) / 2;

        } else {
            s = (maxAlpha - minAlpha) / numberOfLeafs;
            alpha = minAlpha;
        }
        List<ElkNode> successors = RadialUtil.getSuccessors(node);
        if (sorter != null) {
            sorter.initialize(root);
            sorter.sort(successors);
        }
        for (ElkNode child : successors) {
            double numberOfChildLeafs = annulusWedgeCriteria.calculateWedgeSpace(child);
            positionNodes(child, currentRadius + radius, alpha, alpha + s * numberOfChildLeafs, optimalOffset);
            alpha += s * numberOfChildLeafs;
        }
    }

    @Override
    public LayoutProcessorConfiguration<RadialLayoutPhases, ElkNode> getLayoutProcessorConfiguration(
            final ElkNode graph) {
        return null;
    }

}