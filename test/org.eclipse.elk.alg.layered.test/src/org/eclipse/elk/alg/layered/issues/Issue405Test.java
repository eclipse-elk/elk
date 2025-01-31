/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.issues;

import static org.junit.Assert.fail;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.ElkShape;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;
import com.google.common.math.DoubleMath;

/**
 * Test for issue 405.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue405Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        // Note that there are further test cases for #405. However, they cannot be tested in this test case
        // due to assumptions made here regarding the symmetry of computed port and port label positions.
        // They are checked in 'Issue405_2Test'.
        return Lists.newArrayList(
                new ModelResourcePath("tickets/layered/405_differentPortLabelPositionsNSWEInside.elkt"),
                new ModelResourcePath("tickets/layered/405_differentPortLabelPositionsNSWEOutsideNextToPort.elkt"),
                new ModelResourcePath("tickets/layered/405_differentPortLabelPositionsTwoWestOutside.elkt"));
    }

    private static <T extends ElkGraphElement> T getElementWithId(final Collection<T> container, final String id) {
        return container.stream().filter(e -> e.getIdentifier().equals(id)).findFirst().get();
    }

    // Public to be accessible by 'Issue405_2Test'
    public static Pair<List<ElkPort>, List<ElkLabel>> getRelevantPortsAndPortLabels(final ElkNode graph,
            final String rootNodeName) {
        final ElkNode root = getElementWithId(graph.getChildren(), rootNodeName);
        final ElkNode n1 = getElementWithId(root.getChildren(), "N1");
        return Pair.of(n1.getPorts(),
                n1.getPorts().stream().flatMap(p -> p.getLabels().stream()).collect(Collectors.toList()));
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    private static final double DOUBLE_EPSILON = 10e-4;
    
    private <T extends ElkShape> boolean checkPositionEquals(T first, T second) {
        return DoubleMath.fuzzyEquals(first.getX(), second.getX(), DOUBLE_EPSILON)
                && DoubleMath.fuzzyEquals(first.getY(), second.getY(), DOUBLE_EPSILON);
    }
    
    @Test
    public void testPortAndLabelPositionsEqual(final ElkNode graph) {
        // There are four (non-hierarchical) nodes within the diagram for which different layout directions are used.
        // Each node has a number of ports, each of which has a single label.
        // While we do not know in which "order" the ports are placed, we expect the collective of port positions 
        // (and port label positions) to match each other for all of the various layout directions. 
        final List<String> topNodes = Lists.newArrayList("Right", "Down", "Left", "Up");

        // Now check for each subsequent pair of these "layout-direction nodes" that their ports and labels match.
        topNodes.stream()
            .map(name -> getRelevantPortsAndPortLabels(graph, name))
            .reduce((previous, current) -> {
                for (ElkPort port : previous.getFirst()) {
                    // For the previous node's port, find all of the current node's ports that match its position
                    List<ElkPort> matchingPorts = current.getFirst().stream()
                            .filter(other -> checkPositionEquals(port, other))
                            .collect(Collectors.toList());
                    if (matchingPorts.size() != 1) {
                        // We expected there to be exactly one such port!
                        fail(positionsDontMatchErrorMessage(port,  current.getFirst()));
                    }
                }
    
                for (ElkLabel label : previous.getSecond()) {
                    // For the previous node's label, find all of the current node's ports that match its position
                    List<ElkLabel> matchingLabels = current.getSecond().stream()
                            .filter(other -> checkPositionEquals(label, other))
                            .collect(Collectors.toList());
                    if (matchingLabels.size() != 1) {
                        // We expected there to be exactly one such label!
                        fail(positionsDontMatchErrorMessage(label,  current.getSecond()));
                    }
                }
                return current;
            });
    }
    
    private String positionsDontMatchErrorMessage(final ElkGraphElement fixedElement,
            final List<? extends ElkGraphElement> otherElements) {
        
        String message = fixedElement + " of graph " + ElkGraphUtil.containingGraph(fixedElement)
                + " has not exactly one matching element in graph "
                + ElkGraphUtil.containingGraph(otherElements.get(0));
        
        // List all other elements
        for (ElkGraphElement matchingElement : otherElements) {
            message += "\n- " + matchingElement;
        }
        
        return message;
    }

}
