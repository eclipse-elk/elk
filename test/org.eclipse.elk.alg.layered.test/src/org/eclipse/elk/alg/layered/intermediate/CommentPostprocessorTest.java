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
package org.eclipse.elk.alg.layered.intermediate;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.test.framework.LayoutTestRunner;
import org.eclipse.elk.alg.test.framework.annotations.Algorithm;
import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.alg.test.framework.annotations.GraphResourceProvider;
import org.eclipse.elk.alg.test.framework.io.AbstractResourcePath;
import org.eclipse.elk.alg.test.framework.io.ModelResourcePath;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;
import org.junit.runner.RunWith;

import com.google.common.collect.Lists;

/**
 * Basic tests for the {@link CommentPostprocessor}.
 */
@RunWith(LayoutTestRunner.class)
@Algorithm(LayeredOptions.ALGORITHM_ID)
@DefaultConfiguration
public class CommentPostprocessorTest {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sources

    @GraphResourceProvider
    public List<AbstractResourcePath> testGraphs() {
        return Lists.newArrayList(new ModelResourcePath("tests/layered/comment_boxes/**"));
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tests
    
    private static final double TOLERANCE = 0.1;
    
    @Test
    public void testSpacings(final ElkNode graph) {
        // Run our test on all non-comment nodes
        for (ElkNode node : graph.getChildren()) {
            // Ignore comment boxes
            if (node.getProperty(LayeredOptions.COMMENT_BOX)) {
                continue;
            }
            
            // Retrieve applicable spacings
            double commentCommentSpacing = graph.getProperty(LayeredOptions.SPACING_COMMENT_COMMENT);
            double commentNodeSpacing = graph.getProperty(LayeredOptions.SPACING_COMMENT_NODE);
            
            // We go through the top and bottom comments from left to right and keep track where we expect the
            // next comment's x coordinate to be (separately for top and bottom comments)
            double nextTopX = Double.NaN;
            double nextBottomX = Double.NaN;
            
            // We'll also keep track of the width required to place top and bottom comments to check whether they are
            // centred above or below the node
            double topCommentsWidth = 0;
            double bottomCommentsWidth = 0;
            
            for (ElkNode comment : gatherConnectedCommentsSortedByX(node)) {
                if (comment.getY() < node.getY()) {
                    // Top comment; check y coordinate
                    assertEquals("Comment node spacing " + commentNodeSpacing + " violated by comment " + comment,
                            node.getY() - commentNodeSpacing - comment.getHeight(),
                            comment.getY(),
                            TOLERANCE);
                    
                    // Check x coordinate
                    if (!Double.isNaN(nextTopX)) {
                        assertEquals(
                                "Comment comment spacing " + commentCommentSpacing + " violated by comment " + comment,
                                nextTopX,
                                comment.getX(),
                                TOLERANCE);
                    }
                    
                    nextTopX = comment.getX() + comment.getWidth() + commentCommentSpacing;
                    topCommentsWidth += comment.getWidth() + commentCommentSpacing;
                    
                } else {
                    // Bottom comment; check y coordinate
                    assertEquals("Comment node spacing " + commentNodeSpacing + " violated by comment " + comment,
                            node.getY() + node.getHeight() + commentNodeSpacing,
                            comment.getY(),
                            TOLERANCE);
                    
                    // Check x coordinate
                    if (!Double.isNaN(nextBottomX)) {
                        assertEquals(
                                "Comment comment spacing " + commentCommentSpacing + " violated by comment " + comment,
                                nextBottomX,
                                comment.getX(),
                                TOLERANCE);
                    }
                    
                    nextBottomX = comment.getX() + comment.getWidth() + commentCommentSpacing;
                    bottomCommentsWidth += comment.getWidth() + commentCommentSpacing;
                }
            }
            
            // Check if comments are centred
            if (topCommentsWidth > 0) {
                topCommentsWidth -= commentCommentSpacing;
                assertEquals("Top comments not centred above node " + node,
                        node.getX() + node.getWidth() / 2,
                        nextTopX - commentCommentSpacing - (topCommentsWidth / 2),
                        TOLERANCE);
            }
            
            if (bottomCommentsWidth > 0) {
                bottomCommentsWidth -= commentCommentSpacing;
                assertEquals("Bottom comments not centred below node " + node,
                        node.getX() + node.getWidth() / 2,
                        nextBottomX - commentCommentSpacing - (bottomCommentsWidth / 2),
                        TOLERANCE);
            }
        }
    }

    /**
     * Returns the comments connected to the given node, sorted by x coordinate.
     */
    private List<ElkNode> gatherConnectedCommentsSortedByX(final ElkNode node) {
        List<ElkNode> comments = new ArrayList<>();
        
        node.getIncomingEdges().stream()
            .map(inEdge -> inEdge.getSources().get(0))
            .map(shape -> ElkGraphUtil.connectableShapeToNode(shape))
            .filter(comment -> comment.getProperty(LayeredOptions.COMMENT_BOX))
            .forEach(comment -> comments.add(comment));
        
        node.getOutgoingEdges().stream()
            .map(outEdge -> outEdge.getTargets().get(0))
            .map(shape -> ElkGraphUtil.connectableShapeToNode(shape))
            .filter(comment -> comment.getProperty(LayeredOptions.COMMENT_BOX))
            .forEach(comment -> comments.add(comment));
        
        comments.sort((n1, n2) -> Double.compare(n1.getX(), n2.getX()));
        
        return comments;
    }
    
}
