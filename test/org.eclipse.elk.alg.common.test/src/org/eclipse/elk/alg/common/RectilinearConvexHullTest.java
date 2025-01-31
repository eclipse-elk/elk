/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.common;

import java.awt.Color;
import java.awt.Polygon;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.junit.Assert;
import org.junit.Test;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

/**
 * Tests the functionality of the {@link RectilinearConvexHull} class.
 */
public class RectilinearConvexHullTest {

    /**
     * Test the following point set. For each corner of the cross a point to the point set added. 
     * The hull is expected to consist of exactly 3 points per staircase, thus 12 points over all.
     *   
     * <pre>
     *      __
     *   __|  |__
     *  |__    __|
     *     |__|
     * </pre>    
     *     
     */
    @Test
    public void testCross() {
        List<Point> points = Lists.newArrayList();
        points.addAll(
                Arrays.asList(
                        new Point(0d, 1d), new Point(1d, 1d), new Point(1d, 0d), 
                        new Point(2d, 0d), new Point(2d, 1d), new Point(3d, 1d), 
                        new Point(3d, 2d), new Point(2d, 2d), new Point(2d, 3d), 
                        new Point(1d, 3d), new Point(1, 2d), new Point(0d, 2d)));
        
        RectilinearConvexHull rch = RectilinearConvexHull.of(points);
        
        // test that the hull consists of the specified points in the exact order
        Assert.assertEquals(points, rch.getHull());
        
        List<ElkRectangle> expectedRects = Arrays.asList(
                        new ElkRectangle(0, 1, 1, 1), 
                        new ElkRectangle(1, 0, 1, 3), 
                        new ElkRectangle(2, 1, 1, 1));
        List<ElkRectangle> actualRects = rch.splitIntoRectangles();
        
        // we don't care about order here
        Assert.assertTrue(Iterables.elementsEqual(actualRects, expectedRects));
    }
    
    
    /**
     * Test the following point set. 
     *   
     * <pre>
     *     __|__
     *  __|     |__
     *    |_____|
     *       |
     * </pre>    
     *     
     */
    @Test
    public void testLittleRobot() {
        List<Point> points = Lists.newArrayList();
        points.addAll(
                Arrays.asList(
                        new Point(0d, 2d), new Point(1d, 2d), new Point(1d, 1d), new Point(2d, 1d), new Point(2d, 0d), 
                        new Point(2d, 0d), new Point(2d, 1d), new Point(3d, 1d), new Point(3d, 2d), new Point(4d, 2d), 
                        new Point(4d, 2d), new Point(3d, 2d), new Point(3d, 3d), new Point(2d, 3d), new Point(2d, 4d), 
                        new Point(2d, 4d), new Point(2d, 3d), new Point(1d, 3d), new Point(1d, 2d), new Point(0d, 2d)));

        // test hull
        RectilinearConvexHull rch = RectilinearConvexHull.of(points);
        Assert.assertEquals(points, rch.getHull());

        // and rects
        List<ElkRectangle> actualRects = rch.splitIntoRectangles();
        List<ElkRectangle> expectedRects = Arrays.asList(
                new ElkRectangle(0, 2, 1, 0), 
                new ElkRectangle(1, 1, 1, 2), 
                new ElkRectangle(2, 0, 0, 4),
                new ElkRectangle(2, 1, 1, 2),
                new ElkRectangle(3, 2, 1, 0));
        
        Assert.assertTrue(Iterables.elementsEqual(actualRects, expectedRects));
    }
    
    
    /* ----------------
     * Simple JFrame that draws the hull and rects.
     */
    
    public static void main(String[] args) {
        
        List<Point> points = Lists.newArrayList();
        points.addAll(Arrays.asList(new Point(1d, 6d), new Point(3d, 5d), new Point(5d, 2d),
                new Point(7d, 2d), new Point(9d, 4d), new Point(11d, 4d), new Point(11d, 8d),
                new Point(8d, 8d), new Point(8d, 12d), new Point(6d, 12d), new Point(6d, 10d),
                new Point(4d, 6d)));

        Random r = new Random(); 
        for (int i = 0; i < 100; ++i) {
            points.add(new Point((int) (r.nextDouble() * 100), (int) (r.nextDouble() * 100)));
        }

        RectilinearConvexHull rch = RectilinearConvexHull.of(points);
        System.out.println("HULL: " + rch.getHull());

        List<ElkRectangle> rects = rch.splitIntoRectangles();
        System.out.println("Rectangles: " + rects);

        visualize(points, rch);
    }
 
    private static void visualize(final List<Point> points, final RectilinearConvexHull rch) {
        KVector offset = new KVector(10, 10);
        double scale = 10;
        
        
        @SuppressWarnings("serial")
        JPanel p = new JPanel() {

            public void paint(java.awt.Graphics g) {
                super.paint(g);

                g.setColor(Color.DARK_GRAY);
                g.fillPolygon(toPolygon(rch, offset, scale));
                g.setColor(Color.GRAY);
                for (Point p : points) {
                    g.fillRect((int) (p.x * scale + offset.x) - 1,
                            (int) (p.y * scale + offset.y) - 1, 3, 3);
                }
                
                
                g.setColor(Color.PINK);
                for (ElkRectangle r : rch.splitIntoRectangles()) {
                    g.drawRect((int) (r.x * scale + offset.x), (int) ( (r.y) * scale + offset.y),
                            (int) (r.width * scale), (int) (r.height * scale));
                }
                
                g.setColor(Color.GRAY);
            };
        };

        JFrame f = new JFrame();
        f.setSize(1500, 1500);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setVisible(true);
        f.getContentPane().add(p);
        
        try {
            Thread.sleep(10000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
    
    public static Polygon toPolygon(RectilinearConvexHull rch, final KVector offset, final double scale) {
        Polygon poly = new Polygon();
        for (Point pt : rch.getHull()) {
            poly.addPoint((int) (pt.x * scale + offset.x), (int) (pt.y * scale + offset.y));
        }
        return poly;
    }
    
}