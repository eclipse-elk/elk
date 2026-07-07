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
package org.eclipse.elk.alg.disco.debug.views;

import org.eclipse.swt.SWT;
import org.eclipse.swt.graphics.Color;
import org.eclipse.swt.graphics.Font;
import org.eclipse.swt.widgets.Display;

/**
 * Passed to the {@link DisCoGraphRenderer} to configure what the rendering will
 * look like. Clients may subclass this class to override the
 * {@link #initialize(double)} method and create custom instances of the
 * different colors and fonts.
 * 
 * This class started as a copy of GraphRenderingConfigurator .java, commit
 * 71bb8c2f542, 2016-05-07. Changes are made by mic and commented accordingly.
 */
public class DisCoGraphRenderingConfigurator {

    private final Color black;

    /** default font size for nodes. */
    protected static final int NODE_FONT_SIZE = 9;
    /** default font size for ports. */
    protected static final int PORT_FONT_SIZE = 6;
    /** default font size for edges. */
    protected static final int EDGE_FONT_SIZE = 8;
    /** the minimal font height for displaying labels. */
    protected static final int MIN_FONT_HEIGHT = 3;

    // CHECKSTYLEOFF Visibility
    // We think it alright for subclasses to directly set these.

    /** border color for nodes. */
    protected Color nodeBorderColor;
    /** fill color for nodes. */
    protected Color nodeFillColor;
    /** font used for node labels. */
    protected Font nodeLabelFont;
    /** border color used for ports. */
    protected Color portBorderColor;
    /** fill color used for ports. */
    protected Color portFillColor;
    /** font used for port labels. */
    protected Font portLabelFont;
    /** color used for edges. */
    protected Color edgeColor;
    /** font used for edge labels. */
    protected Font edgeLabelFont;
    /** border color for labels. */
    protected Color labelBorderColor;
    /** fill color for labels. */
    protected Color labelFillColor;
    /** text color for labels. */
    protected Color labelTextColor;
    /** by mic: fill color for DCElements. */
    protected Color dcElementFillColor;
    /**
     * 
     */
    protected Color dcElementExternalColor;
    /**
     * 
     */
    protected Color dcElementExternalBorderTextColor;
    /** by mic: border and text color for DCElements. */
    protected Color dcElementBorderTextColor;
    /** by mic: fill color for Polyominoes. */
    protected Color dcPolyominoFillColor;
    /**
     * 
     */
    protected Color dcPolyominoCenterColor;
    /**
     * 
     */
    protected Color dcPolyominoWeaklyBlockedFillColor;
    /**
     * 
     */
    protected Color dcPolyominoWeaklyBlockedBorderTextColor;
    /** by mic: border and text color for Polyominoes. */
    protected Color dcPolyominoBorderTextColor;

    // CHECKSTYLEON Visibility

    /** The display used to create resources. */
    private Display display;

    /**
     * Creates a new instance that uses the given display to create resources.
     * 
     * @param display
     *            the display to create rendering resources on.
     */
    public DisCoGraphRenderingConfigurator(final Display display) {
        this.display = display;
        black = new Color(display, 0, 0, 0);
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
     * Called once by the constructor. Override to create custom resources. Be sure
     * to call the super class implementation unless the configuration methods only
     * return custom resources.
     * 
     * @param scale
     *            Scaling used by the graph renderer. Can be used to adapt font
     *            sizes.
     */
    public void initialize(final double scale) {
        // CHECKSTYLEOFF MagicNumber
        // No sense in introducing constants for the RGB values below...

        // by mic: new color for the kgraph, kind of cyanish
        final Color dcKGraphColor = new Color(display, 0, 115, 127);

        int nodeFontSize = Math.max((int) Math.round(NODE_FONT_SIZE * scale), 2);
        nodeLabelFont = new Font(display, "sans", nodeFontSize, SWT.NORMAL);
        // changed by mic
        nodeBorderColor = dcKGraphColor;
        nodeFillColor = new Color(display, 0, 204, 226);

        // by mic: introducing new colors for DCGraph ans polyominoes.
        dcElementFillColor = new Color(display, 178, 159, 18);
        dcElementExternalColor = new Color(display, 179, 80, 18);
        dcElementBorderTextColor = new Color(display, 127, 113, 13);
        dcElementExternalBorderTextColor = new Color(display, 127, 57, 13);
        dcPolyominoFillColor = new Color(display, 178, 9, 101);
        // dcPolyominoCenterColor = new Color(display, 9, 178, 85);
        dcPolyominoCenterColor = black;
        dcPolyominoWeaklyBlockedFillColor = new Color(display, 170, 9, 178);
        dcPolyominoBorderTextColor = new Color(display, 127, 6, 72);
        dcPolyominoWeaklyBlockedBorderTextColor = new Color(display, 121, 6, 128);

        int portFontSize = Math.max((int) Math.round(PORT_FONT_SIZE * scale), 2);
        portLabelFont = new Font(display, "sans", portFontSize, SWT.NORMAL);
        // changed by mic
        portBorderColor = dcKGraphColor;
        portFillColor = dcKGraphColor;

        int edgeFontSize = Math.max((int) Math.round(EDGE_FONT_SIZE * scale), 2);
        edgeLabelFont = new Font(display, "sans", edgeFontSize, SWT.NORMAL);
        // changed by mic
        edgeColor = dcKGraphColor;

        // changed by mic
        labelBorderColor = dcKGraphColor;
        labelFillColor = null;
        // changed by mic
        labelTextColor = dcKGraphColor;

        // CHECKSTYLEON MagicNumber
    }

    /**
     * Called once by the {@link DisCoGraphRenderer} after the configurator was
     * used.
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

        // by mic: added disposal of newly introduced colors

        if (dcElementFillColor != null) {
            dcElementFillColor.dispose();
        }
        if (dcElementBorderTextColor != null) {
            dcElementBorderTextColor.dispose();
        }
        if (dcPolyominoFillColor != null) {
            dcPolyominoFillColor.dispose();
        }
        if (dcPolyominoBorderTextColor != null) {
            dcPolyominoBorderTextColor.dispose();
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

    // by mic: some new getters for the newly introduced colors.

    public final Color getDCElementFillColor() {
        return dcElementFillColor;
    }

    public final Color getDCElementBorderTextColor() {
        return dcElementBorderTextColor;
    }

    public final Color getDCElementExternalFillColor() {
        return dcElementExternalColor;
    }

    public final Color getDCElementExternalBorderTextColor() {
        return dcElementExternalBorderTextColor;
    }

    public final Color getPolyominoFillColor() {
        return dcPolyominoFillColor;
    }

    public final Color getPolyominoCenterColor() {
        return dcPolyominoCenterColor;
    }

    public final Color getPolyominoBorderTextColor() {
        return dcPolyominoBorderTextColor;
    }

    public final Color getPolyominoWeaklyBlockedColor() {
        return dcPolyominoWeaklyBlockedFillColor;
    }

    public final Color getPolyominoWeaklyBlockedBorderTextColor() {
        return dcPolyominoWeaklyBlockedBorderTextColor;
    }

    /**
     * @return the black
     */
    public final Color getBlack() {
        return black;
    }

}
