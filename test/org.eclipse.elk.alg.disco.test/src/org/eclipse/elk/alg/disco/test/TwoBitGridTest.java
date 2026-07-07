/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.disco.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.eclipse.elk.alg.common.polyomino.structures.TwoBitGrid;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/**
 * Unit test for TwoBitGrid.
 */
public class TwoBitGridTest {
    private TwoBitGrid testGrid;

    /**
     * Some set up code.
     */
    @Before
    public void setUp() {
        final int three = 3;
        final int thirtyfive = 35;
        testGrid = new TwoBitGrid(thirtyfive, three);
    }


    /**
     * Nothing to tear down.
     */
    @After
    public void tearDown() {
    }

    /**
     * Tests whether the dimensions of a grid are set correctly.
     */
    @Test
    public void dimensionTest() {
        final int three = 3;
        final int thirtyfive = 35;
        final int sixtyfour = 64;
        assertEquals(testGrid.getWidth(), thirtyfive);
        assertEquals(testGrid.getHeight(), three);
        TwoBitGrid secondGrid = new TwoBitGrid(sixtyfour, sixtyfour);
        assertEquals(secondGrid.getWidth(), sixtyfour);
        assertEquals(secondGrid.getHeight(), sixtyfour);
    }

    /**
     * Tests whether the predicate {@code inBounds} works.
     */
    @Test
    public void inBoundsTest() {
        final int three = 3;
        final int seventeen = 17;
        final int thirtyfour = 34;
        final int thirtyfive = 35;
        assertTrue(testGrid.inBounds(0, 0));
        assertFalse(testGrid.inBounds(-1, 0));
        assertFalse(testGrid.inBounds(0, -1));
        assertFalse(testGrid.inBounds(-1, -1));
        assertTrue(testGrid.inBounds(thirtyfour, 2));
        assertFalse(testGrid.inBounds(thirtyfour, three));
        assertFalse(testGrid.inBounds(thirtyfive, 2));
        assertFalse(testGrid.inBounds(thirtyfive, three));
        assertTrue(testGrid.inBounds(seventeen, 1));
    }

    /**
     * Checks whether the writing to and reading from specific test cells is working.
     */
    @Test
    public void cellTest() {
        final int thirtythree = 33;
        final int thirtyfour = 34;
        final int sixtyfour = 64;
        
        assertTrue(testGrid.isEmpty(0, 0));
        assertTrue(testGrid.isEmpty(thirtyfour, 2));

        testGrid.setBlocked(0, 0);
        assertTrue(testGrid.isBlocked(0, 0));
        assertFalse(testGrid.isEmpty(0, 0));
        assertFalse(testGrid.isWeaklyBlocked(0, 0));

        testGrid.setWeaklyBlocked(1, 0);
        assertTrue(testGrid.isBlocked(0, 0));
        assertFalse(testGrid.isEmpty(0, 0));
        assertFalse(testGrid.isWeaklyBlocked(0, 0));
        assertFalse(testGrid.isBlocked(1, 0));
        assertFalse(testGrid.isEmpty(1, 0));
        assertTrue(testGrid.isWeaklyBlocked(1, 0));

        testGrid.setEmpty(0, 0);
        assertFalse(testGrid.isBlocked(0, 0));
        assertTrue(testGrid.isEmpty(0, 0));
        assertFalse(testGrid.isWeaklyBlocked(0, 0));
        assertFalse(testGrid.isBlocked(1, 0));
        assertFalse(testGrid.isEmpty(1, 0));
        assertTrue(testGrid.isWeaklyBlocked(1, 0));

        assertTrue(testGrid.isEmpty(thirtythree, 2));

        testGrid.setBlocked(thirtythree, 2);
        assertFalse(testGrid.isBlocked(0, 0));
        assertTrue(testGrid.isEmpty(0, 0));
        assertFalse(testGrid.isWeaklyBlocked(0, 0));
        assertFalse(testGrid.isBlocked(1, 0));
        assertFalse(testGrid.isEmpty(1, 0));
        assertTrue(testGrid.isWeaklyBlocked(1, 0));
        assertFalse(testGrid.isBlocked(1, 0));
        assertFalse(testGrid.isEmpty(1, 0));
        assertTrue(testGrid.isWeaklyBlocked(1, 0));
        assertTrue(testGrid.isBlocked(thirtythree, 2));
        assertFalse(testGrid.isEmpty(thirtythree, 2));
        assertFalse(testGrid.isWeaklyBlocked(thirtythree, 2));

        TwoBitGrid secondGrid = new TwoBitGrid(sixtyfour, sixtyfour);
        for (int yi = 0; yi < secondGrid.getHeight(); yi++) {
            for (int xi = 0; xi < secondGrid.getWidth(); xi++) {
                secondGrid.setBlocked(xi, yi);
            }
        }
        for (int yi = 0; yi < secondGrid.getHeight(); yi++) {
            for (int xi = 0; xi < secondGrid.getWidth(); xi++) {
                assertTrue(secondGrid.isBlocked(xi, yi));
            }
        }

        for (int yi = 0; yi < secondGrid.getHeight(); yi++) {
            for (int xi = 0; xi < secondGrid.getWidth(); xi++) {
                secondGrid.setWeaklyBlocked(xi, yi);
            }
        }
        for (int yi = 0; yi < secondGrid.getHeight(); yi++) {
            for (int xi = 0; xi < secondGrid.getWidth(); xi++) {
                assertTrue(secondGrid.isWeaklyBlocked(xi, yi));
            }
        }

        for (int yi = 0; yi < secondGrid.getHeight(); yi++) {
            for (int xi = 0; xi < secondGrid.getWidth(); xi++) {
                secondGrid.setEmpty(xi, yi);
            }
        }
        for (int yi = 0; yi < secondGrid.getHeight(); yi++) {
            for (int xi = 0; xi < secondGrid.getWidth(); xi++) {
                assertTrue(secondGrid.isEmpty(xi, yi));
            }
        }
    }

}
