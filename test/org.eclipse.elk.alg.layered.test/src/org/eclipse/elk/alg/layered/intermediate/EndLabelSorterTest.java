/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.LabelCell;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.annotations.TestAfterProcessor;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Tests whether the {@link LabelSorter} correctly sorts labels in our test graphs. This test heavily relies on
 * the fact that the text of the labels in the test graphs reflects the expected order, that is, if they are
 * ordered correctly, that order should be equal to a lexicographical order.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration()
public class EndLabelSorterTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tests/layered/labels_ordering/**/"));
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    @TestAfterProcessor(EndLabelSorter.class)
    public void testCorrectOrder(final Object graph) {
        LGraph lGraph = (LGraph) graph;
        
        lGraph.getLayers().stream()
            .flatMap(layer -> layer.getNodes().stream())
            .forEach(node -> testNode(node));
    }
    
    private void testNode(final LNode node) {
        switch (node.getType()) {
        case NORMAL:
            // Build a list of label strings that we can check
            if (node.hasProperty(InternalProperties.END_LABELS)) {
                for (LabelCell labelCell : node.getProperty(InternalProperties.END_LABELS).values()) {
                    List<String> labelStrings = labelCell.getLabels().stream()
                            .map(label -> label.getText())
                            .collect(Collectors.toList());
                    failIfUnsorted(labelStrings);
                }
            }
            break;
            
        case LABEL:
            List<LLabel> labels = node.getProperty(InternalProperties.REPRESENTED_LABELS);
            assertNotNull(labels);
            
            List<String> labelStrings = labels.stream()
                    .map(label -> label.getText())
                    .collect(Collectors.toList());
            failIfUnsorted(labelStrings);
            
            break;
        }
    }

    private void failIfUnsorted(final List<String> strings) {
        if (strings.isEmpty()) {
            return;
        }
        
        Iterator<String> iter = strings.iterator();
        String prev = iter.next();
        while (iter.hasNext()) {
            String curr = iter.next();
            assertTrue(strings.toString(), prev.compareTo(curr) <= 0);
            prev = curr;
        }
    }

}
