/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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

import java.util.List;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.FileNameFilter;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkNode;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Test for issue 515.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class Issue515Test {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(
                new ModelResourcePath("tickets/layered/**").withFilter(new FileNameFilter("515.+\\.elkt")));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests

    @Test
    public void testNoEdgeOverlaps(final ElkNode graph) {
        
        List<ElkRectangle> allNodesAsBoxes = graph.getChildren().stream()
                .map(n -> new ElkRectangle(n.getX(), n.getY(), n.getWidth(), n.getHeight()))
                .collect(Collectors.toList());

        List<KVectorChain> allEdgePaths = graph.getContainedEdges().stream()
                .flatMap(e -> e.getSections().stream())
                .map(s -> ElkUtil.createVectorChain(s))
                .collect(Collectors.toList());
        
        allEdgePaths.forEach(path -> {
            // For each edge segment (start, end), check if it intersects any box
            path.stream().reduce((start, end) -> {
                allNodesAsBoxes.forEach(box -> {
                    assertTrue("Edge path overlaps with node.", !ElkMath.intersects(box, start, end));
                });
                return end;
            });
        });
    }

}
