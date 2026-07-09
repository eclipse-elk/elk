/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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
package org.eclipse.elk.alg.layered.issues;

import static org.junit.Assert.assertTrue;

import java.awt.geom.Rectangle2D;
import java.util.List;
import java.util.Optional;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 433.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue433Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tickets/layered/433_graphSizeWithSelfLoopLabel.elkt"));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    /**
     * Check that the self loop label in the graph is fully contained in the graph.
     */
    @Test
    public void testProperLabelCoordinates(final ElkNode graph) {
        Rectangle2D.Double graphRect = new Rectangle2D.Double(0, 0, graph.getWidth(), graph.getHeight());
        
        Optional<Rectangle2D.Double> labelRect = graph.getContainedEdges().stream()
                .flatMap(edge -> edge.getLabels().stream())
                .map(label -> new Rectangle2D.Double(label.getX(), label.getY(), label.getWidth(), label.getHeight()))
                .findFirst();
        
        assertTrue("Self loop label is not entirely inside the graph."
                + "\nGraph rect: " + graphRect + "\nLabel rect: " + labelRect.get(),
                graphRect.contains(labelRect.get()));
    }

}
