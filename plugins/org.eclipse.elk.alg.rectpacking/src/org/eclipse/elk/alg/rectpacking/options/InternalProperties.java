/*******************************************************************************
 * Copyright (c) 2022 Kiel Univerity and others.
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
package org.eclipse.elk.alg.rectpacking.options;

import java.util.List;

import org.eclipse.elk.alg.rectpacking.util.RectRow;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;

/**
 * Container for property definitions for internal use of the algorithm. These properties should
 * not be accessed from outside.
 */
public final class InternalProperties {

    /** Additional height of the drawing that is to be added during expansion. */
    public static final IProperty<Double> ADDITIONAL_HEIGHT = new Property<Double>("additionalHeight");
    /** The height of the parent. */
    public static final IProperty<Double> DRAWING_HEIGHT = new Property<Double>("drawingHeight");
    /** The width of the parent. */
    public static final IProperty<Double> DRAWING_WIDTH = new Property<Double>("drawingWidth");

    /** The minimum height of the parent. */
    public static final IProperty<Double> MIN_HEIGHT = new Property<Double>("minHeight");
    /** The minimum width of the parent. */
    public static final IProperty<Double> MIN_WIDTH = new Property<Double>("minWidth");
    
    /** The rows formed by the simple placement. */
    public static final IProperty<List<RectRow>> ROWS = new Property<List<RectRow>>("rows");

    /** The target width of the graph. */
    public static final IProperty<Double> TARGET_WIDTH = new Property<Double>("targetWidth");

    /**
     * The minimum width that has to be added to the target width to allow an element from the next row to be placed
     * the previous row.
     */
    public static final IProperty<Double> MIN_ROW_INCREASE = new Property<Double>("minRowIncrease", 0.0);
    /**
     * The maximum width that has to be added to the target width to allow an element from the next row to be placed
     * the previous row.
     */
    public static final IProperty<Double> MAX_ROW_INCREASE = new Property<Double>("maxRowIncrease", 0.0);

    /**
     * The minimum width that has to be removed from the target width to not allow the last element of a row in it.
     */
    public static final IProperty<Double> MIN_ROW_DECREASE = new Property<Double>("minRowDecrease", 0.0);

    /**
     * The maximum width that has to be removed from the target width to not allow the last element of a row in it.
     */
    public static final IProperty<Double> MAX_ROW_DECREASE = new Property<Double>("maxRowDecrease", 0.0);

    /**
     * hidden default constructor.
     */
    private InternalProperties() {

    }
}
