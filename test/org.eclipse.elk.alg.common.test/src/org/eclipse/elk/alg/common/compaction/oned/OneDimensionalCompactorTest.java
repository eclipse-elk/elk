/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.compaction.oned;

import static org.junit.Assert.assertEquals;

import java.util.EnumSet;
import java.util.Set;

import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.options.Direction;
import org.junit.Test;

/**
 * Several test for the {@link OneDimensionalCompactor}.
 */
public class OneDimensionalCompactorTest {

    public static final double EPSILON = 0.0001d;
    public static final double SPACING = 5;
    
    public static final ISpacingsHandler TEST_SPACING_HANDLER = new ISpacingsHandler() {
        @Override
        public double getVerticalSpacing(CNode cNode1, CNode cNode2) {
            return SPACING;
        }
        @Override
        public double getHorizontalSpacing(CNode cNode1, CNode cNode2) {
            return SPACING;
        }
    }; 
    
    @Test
    public void testLeftCompaction() {
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode left = CNode.of()
                          .hitbox(new ElkRectangle(0, 0, 20, 20))
                          .create(graph);
        CNode right = CNode.of()
                           .hitbox(new ElkRectangle(30, 0, 20, 20))
                           .create(graph);
        
        compacter(graph)
                 .changeDirection(Direction.LEFT)
                 .compact()
                 .finish();
        
        assertEquals(0, left.hitbox.x, EPSILON);
        assertEquals(20, right.hitbox.x, EPSILON);
    }
    
    @Test 
    public void testLeftCompactionEqualYCoordinate() {
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode top = CNode.of()
                .hitbox(new ElkRectangle(0, 0, 20, 20))
                .create(graph);
        CNode bot = CNode.of()
                .hitbox(new ElkRectangle(30, 20, 20, 20))
                .create(graph);
        
        compacter(graph)
                 .setConstraintAlgorithm(OneDimensionalCompactor.SCANLINE_CONSTRAINTS)
                 .changeDirection(Direction.LEFT)
                 .compact()
                 .finish();
        
        assertEquals(0, top.hitbox.x, EPSILON);
        assertEquals(0, bot.hitbox.x, EPSILON);
    }
    
    @Test
    public void testLeftCompactionSpacingAware() {
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode left = CNode.of()
                          .hitbox(new ElkRectangle(0, 0, 20, 20))
                          .create(graph);
        CNode right = CNode.of()
                           .hitbox(new ElkRectangle(30, 20 + SPACING - 1, 20, 20))
                           .create(graph);
        
        compacter(graph)
                 .setSpacingsHandler(TEST_SPACING_HANDLER)
                 .changeDirection(Direction.LEFT)
                 .compact()
                 .finish();
        
        assertEquals(0, left.hitbox.x, EPSILON);
        assertEquals(25, right.hitbox.x, EPSILON);
    }
    
    @Test
    public void testLeftCompactionSpacingAware2() {
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode left = CNode.of()
                          .hitbox(new ElkRectangle(0, 0, 20, 20))
                          .create(graph);
        CNode right = CNode.of()
                           .hitbox(new ElkRectangle(30, 20 + SPACING + 1, 20, 20))
                           .create(graph);
        
        compacter(graph)
                 .setSpacingsHandler(TEST_SPACING_HANDLER)
                 .changeDirection(Direction.LEFT)
                 .compact()
                 .finish();
        
        assertEquals(0, left.hitbox.x, EPSILON);
        assertEquals(0, right.hitbox.x, EPSILON);
    }
    
    
    @Test
    public void testRightCompaction() {
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode left = CNode.of()
                          .hitbox(new ElkRectangle(0, 0, 20, 20))
                          .create(graph);
        CNode right = CNode.of()
                           .hitbox(new ElkRectangle(30, 0, 20, 20))
                           .create(graph);
        
        compacter(graph)
                 .changeDirection(Direction.RIGHT)
                 .compact()
                 .finish();
        
        assertEquals(10, left.hitbox.x, EPSILON);
        assertEquals(30, right.hitbox.x, EPSILON);
    }
    
    @Test
    public void testUpCompaction() {
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode upper = CNode.of()
                .hitbox(new ElkRectangle(0, 0, 20, 20))
                .create(graph);
        CNode lower = CNode.of()
                .hitbox(new ElkRectangle(0, 30, 20, 20))
                .create(graph);
        
        compacter(graph)
                 .changeDirection(Direction.UP)
                 .compact()
                 .finish();
        
        assertEquals(0, upper.hitbox.y, EPSILON);
        assertEquals(20, lower.hitbox.y, EPSILON);
    }
    
    @Test
    public void testDownCompaction() {
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode upper = CNode.of()
                .hitbox(new ElkRectangle(0, 0, 20, 20))
                .create(graph);
        CNode lower = CNode.of()
                .hitbox(new ElkRectangle(0, 30, 20, 20))
                .create(graph);
        
        compacter(graph)
                 .changeDirection(Direction.DOWN)
                 .compact()
                 .finish();
        
        assertEquals(10, upper.hitbox.y, EPSILON);
        assertEquals(30, lower.hitbox.y, EPSILON);
    }
    
    /**
     * The connection below indicates a grouping, not an edge.
     * 
     *   +--+         
     *   |  |     +--+
     *   +--+     |  |
     *            +--+
     *             |  
     *          +--+  
     *          |  |  
     *          +--+ 
     */
    @Test 
    public void testLeftGroupCompaction() {
        
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode left = CNode.of()
                .hitbox(new ElkRectangle(0, 0, 20, 20))
                .create(graph);
        CNode upperRight = CNode.of()
                .hitbox(new ElkRectangle(40, 5, 20, 20))
                .create(graph);
        CNode lowerRight = CNode.of()
                .hitbox(new ElkRectangle(30, 25, 20, 20))
                .create(graph);
        
        CGroup.of()
                .nodes(upperRight, lowerRight)
                .create(graph);
        
        compacter(graph)
                .changeDirection(Direction.LEFT)
                .compact()
                .finish();
        
        assertEquals(0, left.hitbox.x, EPSILON);
        assertEquals(20, upperRight.hitbox.x, EPSILON);
        assertEquals(10, lowerRight.hitbox.x, EPSILON);
    }
    
    /**
     * The connection below indicates a grouping, not an edge.
     * 
     *   +--+         
     *   |  |     +--+
     *   +--+     |  |
     *     |      +--+
     *     |          
     *     +--+  
     *     |  |  
     *     +--+ 
     */
    @Test 
    public void testRightGroupCompaction() {
        
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode left = CNode.of()
                .hitbox(new ElkRectangle(0, 5, 20, 20))
                .create(graph);
        CNode upperRight = CNode.of()
                .hitbox(new ElkRectangle(40, 0, 20, 20))
                .create(graph);
        CNode lowerRight = CNode.of()
                .hitbox(new ElkRectangle(10, 25, 20, 20))
                .create(graph);
        
        CGroup.of()
            .nodes(left, lowerRight)
            .create(graph);
        
        compacter(graph)
                .changeDirection(Direction.RIGHT)
                .compact()
                .finish();
        
        assertEquals(20, left.hitbox.x, EPSILON);
        assertEquals(40, upperRight.hitbox.x, EPSILON);
        assertEquals(30, lowerRight.hitbox.x, EPSILON);
    }
    
    /**
     * The connection below indicates a grouping, not an edge.
     * 
     *   +--+         
     *   |  |     
     *   +--+     +--+
     *            |  |
     *     +--+---+--+ 
     *     |  |  
     *     +--+
     */
    @Test 
    public void testUpGroupCompaction() {
        
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode upperLeft = CNode.of()
                .hitbox(new ElkRectangle(0, 0, 20, 20))
                .create(graph);
        CNode lowerLeft = CNode.of()
                .hitbox(new ElkRectangle(5, 40, 20, 20))
                .create(graph);
        CNode right = CNode.of()
                .hitbox(new ElkRectangle(25, 30, 20, 20))
                .create(graph);
        
        CGroup.of()
            .nodes(lowerLeft, right)
            .create(graph);
        
        compacter(graph)
                .changeDirection(Direction.UP)
                .compact()
                .finish();
        
        assertEquals(0, upperLeft.hitbox.y, EPSILON);
        assertEquals(20, lowerLeft.hitbox.y, EPSILON);
        assertEquals(10, right.hitbox.y, EPSILON);
    }
    
    /**
     * The connection below indicates a grouping, not an edge.
     * 
     *   +--+         
     *   |  |     
     *   +--+-----+--+
     *            |  |
     *     +--+   +--+ 
     *     |  |  
     *     +--+ 
     */
    @Test 
    public void testDownGroupCompaction() {
        
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode upperLeft = CNode.of()
                .hitbox(new ElkRectangle(0, 0, 20, 20))
                .create(graph);
        CNode lowerLeft = CNode.of()
                .hitbox(new ElkRectangle(5, 40, 10, 20))
                .create(graph);
        CNode right = CNode.of()
                .hitbox(new ElkRectangle(25, 10, 20, 20))
                .create(graph);
        
        CGroup.of()
            .nodes(upperLeft, right)
            .create(graph);
        
        compacter(graph)
                .changeDirection(Direction.DOWN)
                .compact()
                .finish();
        
        assertEquals(20, upperLeft.hitbox.y, EPSILON);
        assertEquals(40, lowerLeft.hitbox.y, EPSILON);
        assertEquals(30, right.hitbox.y, EPSILON);
    }
    
    /**
     * Two separate groups, three nodes each.
     * 
     * Care has to be taken when calculating the initial position of a {@link CGroup} not to include
     * any spacing between the group nodes. (Group nodes are rigidly fixed relative to each other
     * and thus no other spacing should be applied in between them). This test asserts that the
     * initial position of the groups is as expected.
     */
    @Test
    public void testNoSpacingAppliedWithinGroups() {
        
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode one = CNode.of()
                .hitbox(new ElkRectangle(0, 0, 20, 20))
                .create(graph);
        CNode two = CNode.of()
                .hitbox(new ElkRectangle(20, 10, 20, 20))
                .create(graph);
        CNode three = CNode.of()
                .hitbox(new ElkRectangle(40, 20, 20, 20))
                .create(graph);
        
        CGroup.of()
            .nodes(one, two, three)
            .create(graph);
        
        CNode four = CNode.of()
                .hitbox(new ElkRectangle(22, 80, 20, 20))
                .create(graph);
        CNode five = CNode.of()
                .hitbox(new ElkRectangle(42, 90, 20, 20))
                .create(graph);
        CNode six = CNode.of()
                .hitbox(new ElkRectangle(62, 100, 20, 20))
                .create(graph);
        
        CGroup.of()
            .nodes(four, five, six)
            .create(graph);
        
        compacter(graph)
            .setSpacingsHandler(TEST_SPACING_HANDLER)
            .changeDirection(Direction.LEFT)
            .compact()
            .changeDirection(Direction.RIGHT)
            .compact()
            .changeDirection(Direction.UP)
            .compact()
            .changeDirection(Direction.DOWN)
            .compact()
            .finish();

        assertEquals(0, one.hitbox.x, EPSILON);
        assertEquals(20, two.hitbox.x, EPSILON);
        assertEquals(40, three.hitbox.x, EPSILON);
        assertEquals(0, four.hitbox.x, EPSILON);
        assertEquals(20, five.hitbox.x, EPSILON);
        assertEquals(40, six.hitbox.x, EPSILON);
        
        assertEquals(0, one.hitbox.y, EPSILON);
        assertEquals(10, two.hitbox.y, EPSILON);
        assertEquals(20, three.hitbox.y, EPSILON);
        assertEquals(35, four.hitbox.y, EPSILON);
        assertEquals(45, five.hitbox.y, EPSILON);
        assertEquals(55, six.hitbox.y, EPSILON);
    }
    
    /* --------------------------------------------------
     * Testing subsequent calls with different directions
     * -------------------------------------------------- */
    
    @Test 
    public void testSubsequentDirectionsCompaction() {
        
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode one = CNode.of()
                .hitbox(new ElkRectangle(0, 0, 20, 20))
                .create(graph);
        CNode two = CNode.of()
                .hitbox(new ElkRectangle(25, 0, 20, 20))
                .create(graph);
        CNode three = CNode.of()
                .hitbox(new ElkRectangle(0, 25, 20, 20))
                .create(graph);
        CNode four = CNode.of()
                .hitbox(new ElkRectangle(25, 25, 20, 20))
                .create(graph);
        
        Set<Direction> directions = EnumSet.of(Direction.LEFT, Direction.RIGHT, Direction.UP, Direction.DOWN);
        
        // subsequently apply all combinations of four subsequent compaction steps
        for (Direction d1 : directions) {
            for (Direction d2 : directions) {
                for (Direction d3 : directions) {
                    for (Direction d4 : directions) {
                        
                        compacter(graph)
                            .setSpacingsHandler(TEST_SPACING_HANDLER)
                            .changeDirection(d1)
                            .compact()
                            .changeDirection(d2)
                            .compact()
                            .changeDirection(d3)
                            .compact()
                            .changeDirection(d4)
                            .compact()
                            .finish();
                
                        // the way we modeled the graph, every node should stay where it is 
                        String currentDirections = d1 + " " + d2 + " " + d3 + " " + d4;
                        assertEquals(currentDirections, 0, one.hitbox.x, EPSILON);
                        assertEquals(currentDirections, 0, one.hitbox.y, EPSILON);
                        
                        assertEquals(currentDirections, 25, two.hitbox.x, EPSILON);
                        assertEquals(currentDirections, 0, two.hitbox.y, EPSILON);
                        
                        assertEquals(currentDirections, 0, three.hitbox.x, EPSILON);
                        assertEquals(currentDirections, 25, three.hitbox.y, EPSILON);
                        
                        assertEquals(currentDirections, 25, four.hitbox.x, EPSILON);
                        assertEquals(currentDirections, 25, four.hitbox.y, EPSILON);
                    }
                }
            }
        }
        
    }
    
    /* --------------------------------------------------
     *        Testing different kinds of spacings.
     * --------------------------------------------------*/
    
    /*
     * We support different spacing values between nodes and edges in vertical and horizontal
     * direction. Thereby "vertical" and "horizontal" depends on the compaction direction.
     * Furthermore, we support individual spacings between any pair of nodes.
     */
    
    /**
     * Horizontal spacing should be preserved when compaction leftwards or rightwards.
     */
    @Test
    public void testHorizontalSpacings() {

        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        // test horizontal direction
        CNode one = CNode.of().hitbox(new ElkRectangle(0, 0, 20, 20)).create(graph);
        CNode two = CNode.of().hitbox(new ElkRectangle(50, 0, 20, 20)).create(graph);
        CNode three = CNode.of().hitbox(new ElkRectangle(150, 0, 20, 20)).create(graph);
        
        ISpacingsHandler spacingHandler = new ISpacingsHandler() {
            @Override
            public double getVerticalSpacing(CNode cNode1, CNode cNode2) { return SPACING; }
            @Override
            public double getHorizontalSpacing(CNode cNode1, CNode cNode2) {
                // order realizes a "max" function
                if (cNode1 == three || cNode2 == three) return 10;
                if (cNode1 == two || cNode2 == two) return 7;
                if (cNode1 == one || cNode2 == one) return 5;
                return SPACING;
            }
        };
        
        compacter(graph)
            .setSpacingsHandler(spacingHandler)
            .changeDirection(Direction.LEFT)
            .compact()
            .finish();
        
        assertEquals(0, one.hitbox.x, EPSILON);
        assertEquals(27, two.hitbox.x, EPSILON);
        assertEquals(57, three.hitbox.x, EPSILON);
        
        compacter(graph)
            .setSpacingsHandler(spacingHandler)
            .changeDirection(Direction.RIGHT)
            .compact()
            .finish();
    
        assertEquals(0, one.hitbox.x, EPSILON);
        assertEquals(27, two.hitbox.x, EPSILON);
        assertEquals(57, three.hitbox.x, EPSILON);

        compacter(graph)
            .setSpacingsHandler(spacingHandler)
            .changeDirection(Direction.LEFT)
            .compact()
            .finish();
    
        assertEquals(0, one.hitbox.x, EPSILON);
        assertEquals(27, two.hitbox.x, EPSILON);
        assertEquals(57, three.hitbox.x, EPSILON);
    }
    
    /**
     * "Two" is "center node". Node "one" is supposed to be blocked by "two", while "three" can move freely.
     */
    @Test
    public void testVerticalSpacingDuringHorizontalCompaction() {
        
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode one = CNode.of().hitbox(new ElkRectangle(150, 11, 20, 20)).create(graph);
        CNode two = CNode.of().hitbox(new ElkRectangle(0, 40, 20, 20)).create(graph);
        CNode three = CNode.of().hitbox(new ElkRectangle(150, 76, 20, 20)).create(graph);
        
        ISpacingsHandler spacingHandler = new ISpacingsHandler() {
            
            @Override
            public double getVerticalSpacing(CNode cNode1, CNode cNode2) {
                // order realizes a "max" function
                if (cNode1 == three || cNode2 == three) return 15;
                if (cNode1 == one || cNode2 == one) return 10;
                if (cNode1 == two || cNode2 == two) return 5;
                return 0;
            }
            
            @Override
            public double getHorizontalSpacing(CNode cNode1, CNode cNode2) { return 0; }
        };
        
        compacter(graph)
            .setSpacingsHandler(spacingHandler)
            .changeDirection(Direction.LEFT)
            .compact()
            .finish();
        
        assertEquals(20, one.hitbox.x, EPSILON);
        assertEquals(0, two.hitbox.x, EPSILON);
        assertEquals(0, three.hitbox.x, EPSILON);
    }
    
    /**
     * Vertical spacing should be preserved when compaction upwards or downwards.
     */
    @Test
    public void testVerticalSpacings() {

        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        // test horizontal direction
        CNode one = CNode.of().hitbox(new ElkRectangle(0, 0, 20, 20)).create(graph);
        CNode two = CNode.of().hitbox(new ElkRectangle(0, 50, 20, 20)).create(graph);
        CNode three = CNode.of().hitbox(new ElkRectangle(0, 150, 20, 20)).create(graph);
        
        ISpacingsHandler spacingHandler = new ISpacingsHandler() {
            
            @Override
            public double getVerticalSpacing(CNode cNode1, CNode cNode2) {
                // order realizes a "max" function
                if (cNode1 == three || cNode2 == three) return 10;
                if (cNode1 == two || cNode2 == two) return 7;
                if (cNode1 == one || cNode2 == one) return 5;
                return 0;
            }
            
            @Override
            public double getHorizontalSpacing(CNode cNode1, CNode cNode2) { return 0; }
        };
        
        compacter(graph)
            .setSpacingsHandler(spacingHandler)
            .changeDirection(Direction.UP)
            .compact()
            .finish();
        
        assertEquals(0, one.hitbox.y, EPSILON);
        assertEquals(27, two.hitbox.y, EPSILON);
        assertEquals(57, three.hitbox.y, EPSILON);
        
        compacter(graph)
            .setSpacingsHandler(spacingHandler)
            .changeDirection(Direction.DOWN)
            .compact()
            .finish();
        
        assertEquals(0, one.hitbox.y, EPSILON);
        assertEquals(27, two.hitbox.y, EPSILON);
        assertEquals(57, three.hitbox.y, EPSILON);
        
        compacter(graph)
            .setSpacingsHandler(spacingHandler)
            .changeDirection(Direction.UP)
            .compact()
            .finish();
        
        assertEquals(0, one.hitbox.y, EPSILON);
        assertEquals(27, two.hitbox.y, EPSILON);
        assertEquals(57, three.hitbox.y, EPSILON);
    }
    
    /**
     * "Two" is "center node". Node "one" is supposed to be blocked by "two", while "three" can move freely.
     */
    @Test
    public void testHorizontalSpacingDuringVerticalCompaction() {
        
        CGraph graph = new CGraph(EnumSet.allOf(Direction.class));
        
        CNode one = CNode.of().hitbox(new ElkRectangle(16, 150, 20, 20)).create(graph);
        CNode two = CNode.of().hitbox(new ElkRectangle(40, 0, 20, 20)).create(graph);
        CNode three = CNode.of().hitbox(new ElkRectangle(76, 150, 20, 20)).create(graph);
        
        ISpacingsHandler spacingHandler = new ISpacingsHandler() {
            
            @Override
            public double getVerticalSpacing(CNode cNode1, CNode cNode2) { return 0; }
            
            @Override
            public double getHorizontalSpacing(CNode cNode1, CNode cNode2) { 
                // order realizes a "max" function
                if (cNode1 == three || cNode2 == three) return 15;
                if (cNode1 == two || cNode2 == two) return 10;
                if (cNode1 == one || cNode2 == one) return 5;
                return 0; 
            }
        };
        
        compacter(graph)
            .setSpacingsHandler(spacingHandler)
            .changeDirection(Direction.UP)
            .compact()
            .finish();
        
        assertEquals(20, one.hitbox.y, EPSILON);
        assertEquals(0, two.hitbox.y, EPSILON);
        assertEquals(0, three.hitbox.y, EPSILON);
    }
    
    
    ////////////////////////////////// Internal Convenience API //////////////////////////////////
    
    private OneDimensionalCompactor compacter(final CGraph graph) {
        return new OneDimensionalCompactor(graph)
            .setConstraintAlgorithm(OneDimensionalCompactor.QUADRATIC_CONSTRAINTS);
        
        // TODO test the other constraint algorithm as well
    }

}
