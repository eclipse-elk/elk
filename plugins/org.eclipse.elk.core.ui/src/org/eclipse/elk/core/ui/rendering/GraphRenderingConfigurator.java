/*******************************************************************************
 * Copyright (c) 2016, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.ui.rendering;

import org.eclipse.swt.SWT;
import org.eclipse.swt.graphics.Color;
import org.eclipse.swt.graphics.Font;
import org.eclipse.swt.widgets.Display;

/**
 * Passed to the {@link GraphRenderer} to configure what the rendering will look like. Clients may
 * subclass this class to override the {@link #initialize(double)} method and create custom instances
 * of the different colors and fonts.
 */
public class GraphRenderingConfigurator {
    
    /** default font size for nodes. */
    protected static final int NODE_FONT_SIZE = 9;
    /** default font size for ports. */
    protected static final int PORT_FONT_SIZE = 6;
    /** default font size for edges. */
    protected static final int EDGE_FONT_SIZE = 8;
    /** the minimal font height for displaying labels. */
    protected static final int MIN_FONT_HEIGHT = 3;
    
    // CHECKSTYLEOFF Visibility
    // We think it is alright for subclasses to directly set these.
    
    /** border color for nodes. */
    protected Color nodeBorderColor;
    /** fill color for nodes. */
    protected Color nodeFillColor;
    /** font used for node labels. */
    protected Font nodeLabelFont;
    private int nodeFontSize;
    /** border color used for ports. */
    protected Color portBorderColor;
    /** fill color used for ports. */
    protected Color portFillColor;
    /** font used for port labels. */
    protected Font portLabelFont;
    private int portFontSize;
    /** color used for edges. */
    protected Color edgeColor;
    /** font used for edge labels. */
    protected Font edgeLabelFont;
    private int edgeFontSize;
    /** border color for labels. */
    protected Color labelBorderColor;
    /** fill color for labels. */
    protected Color labelFillColor;
    /** text color for labels. */
    protected Color labelTextColor;
    /** color for area out of bounds. */
    protected Color outOfBoundsColor;
    /** color for the parent node. */
    protected Color rootNodeColor;
    
    // CHECKSTYLEON Visibility
    
    /** The display used to create resources. */
    private Display display;
    
    
    /**
     * Creates a new instance that uses the given display to create resources.
     * 
     * @param display
     *            the display to create rendering resources on.
     */
    public GraphRenderingConfigurator(final Display display) {
        this.display = display;
        initialize();
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Getters / Setters
    
    /**
     * Returns the display this configurator was created for.
     * 
     * @return the display.
     */
    public final Display getDisplay() {
        return display;
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Resource Management
    
    /**
     * Called once by the constructor. Override to create custom resources. Be sure to call the
     * super class implementation unless the configuration methods only return custom resources.
     */
    protected void initialize() {
        // CHECKSTYLEOFF MagicNumber
        // No sense in introducing constants for the RGB values below...
        
        nodeLabelFont = new Font(display, "sans", NODE_FONT_SIZE, SWT.NORMAL);
        nodeBorderColor = new Color(display, 2, 15, 3);
        nodeFillColor = new Color(display, 168, 220, 190);
        
        portLabelFont = new Font(display, "sans", PORT_FONT_SIZE, SWT.NORMAL);
        portBorderColor = new Color(display, 2, 9, 40);
        portFillColor = new Color(display, 2, 9, 40);
        
        edgeLabelFont = new Font(display, "sans", EDGE_FONT_SIZE, SWT.NORMAL);
        edgeColor = new Color(display, 23, 36, 54);
        
        labelBorderColor = new Color(display, 63, 117, 67);
        labelFillColor = null;
        labelTextColor = new Color(display, 2, 15, 3);
        
        outOfBoundsColor = new Color(display, 255, 205, 210);
        rootNodeColor = new Color(display, 240, 240, 240);
        // CHECKSTYLEON MagicNumber
    }
    
    /**
     * Called once by the {@link GraphRenderer} after the configurator was used.
     */
    public final void dispose() {
        if (nodeBorderColor != null) {
            nodeBorderColor.dispose();
        }
        
        if (nodeFillColor != null) {
            nodeFillColor.dispose();
        }
        
        if (nodeLabelFont != null) {
            nodeLabelFont.dispose();
        }
        
        if (labelBorderColor != null) {
            labelBorderColor.dispose();
        }
        
        if (labelTextColor != null) {
            labelTextColor.dispose();
        }
        
        if (portBorderColor != null) {
            portBorderColor.dispose();
        }
        
        if (portFillColor != null) {
            portFillColor.dispose();
        }
        
        if (portLabelFont != null) {
            portLabelFont.dispose();
        }
        
        if (edgeColor != null) {
            edgeColor.dispose();
        }
        
        if (edgeLabelFont != null) {
            edgeLabelFont.dispose();
        }
        
        if (outOfBoundsColor != null) {
            outOfBoundsColor.dispose();
        }
        
        if (rootNodeColor != null) {
            rootNodeColor.dispose();
        }
    }
    
    /**
     * Change the scale of the graph rendering.
     * 
     * @param scale the scale value for font sizes
     */
    public void setScale(final double scale) {
        int newNodeFontSize = Math.max((int) Math.round(NODE_FONT_SIZE * scale), 2);
        if (newNodeFontSize != nodeFontSize) {
            if (nodeLabelFont != null) {
                nodeLabelFont.dispose();
            }
            nodeLabelFont = new Font(display, "sans", newNodeFontSize, SWT.NORMAL);
            nodeFontSize = newNodeFontSize;
        }
        
        int newPortFontSize = Math.max((int) Math.round(PORT_FONT_SIZE * scale), 2);
        if (newPortFontSize != portFontSize) {
            if (portLabelFont != null) {
                portLabelFont.dispose();
            }
            portLabelFont = new Font(display, "sans", newPortFontSize, SWT.NORMAL);
            portFontSize = newPortFontSize;
        }
        
        int newEdgeFontSize = Math.max((int) Math.round(EDGE_FONT_SIZE * scale), 2);
        if (newEdgeFontSize != edgeFontSize) {
            if (edgeLabelFont != null) {
                edgeLabelFont.dispose();
            }
            edgeLabelFont = new Font(display, "sans", newEdgeFontSize, SWT.NORMAL);
            edgeFontSize = newEdgeFontSize;
        }
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Getters
    
    // CHECKSTYLEOFF Javadoc
    // Just a simple set of sufficiently well named getters...

    public final Color getNodeBorderColor() {
        return nodeBorderColor;
    }

    public final Color getNodeFillColor() {
        return nodeFillColor;
    }

    public final Font getNodeLabelFont() {
        return nodeLabelFont;
    }

    public final Color getPortBorderColor() {
        return portBorderColor;
    }

    public final Color getPortFillColor() {
        return portFillColor;
    }

    public final Font getPortLabelFont() {
        return portLabelFont;
    }

    public final Color getEdgeColor() {
        return edgeColor;
    }

    public final Font getEdgeLabelFont() {
        return edgeLabelFont;
    }

    public final Color getLabelBorderColor() {
        return labelBorderColor;
    }

    public final Color getLabelFillColor() {
        return labelFillColor;
    }

    public final Color getLabelTextColor() {
        return labelTextColor;
    }
    
    public final Color getOutOfBoundsColor() {
        return outOfBoundsColor;
    }
    
    public final Color getRootNodeColor() {
        return rootNodeColor;
    }
}
