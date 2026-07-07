/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.disco.debug.views;

import org.eclipse.elk.graph.ElkNode;
import org.eclipse.swt.SWT;
import org.eclipse.swt.events.PaintEvent;
import org.eclipse.swt.events.PaintListener;
import org.eclipse.swt.graphics.Color;
import org.eclipse.swt.graphics.Point;
import org.eclipse.swt.graphics.Rectangle;
import org.eclipse.swt.widgets.Canvas;
import org.eclipse.swt.widgets.Composite;

/**
 * A canvas that is able to paint ELK layout graphs, DCGraphs and polyominoes.
 * Colors and fonts used for painting can be customized by supplying a subclass
 * of {@link DisCoGraphRenderingConfigurator}.
 * 
 * This class started as a copy of GraphRenderingCanvas.java, commit
 * 71bb8c2f542, 2016-05-07. Changes are made by mic and commented accordingly.
 */
public class DisCoGraphRenderingCanvas extends Canvas implements PaintListener {

    /** the painted layout graph. */
    private ElkNode layoutGraph;
    /** the graph renderer used for painting. Replaced by mic */
    private DisCoGraphRenderer graphRenderer;

    /** background color. */
    private Color backgroundColor;

    /**
     * Creates a layout graph canvas using the default drawing style.
     * 
     * @param parent
     *            the parent widget
     */
    public DisCoGraphRenderingCanvas(final Composite parent) {
        super(parent, SWT.NONE);
        addPaintListener(this);
        graphRenderer = new DisCoGraphRenderer(new DisCoGraphRenderingConfigurator(parent.getDisplay()));

        // SUPPRESS CHECKSTYLE NEXT MagicNumber
        backgroundColor = new Color(parent.getDisplay(), 255, 255, 255);
        setBackground(backgroundColor);
    }

    /**
     * Creates a layout graph canvas.
     * 
     * @param parent
     *            the parent widget
     * @param configurator
     *            the rendering configurator that defines the colors and fonts
     *            used for drawing. The configurator must have been initialized
     *            with the same display this component is created for.
     * @throws IllegalArgumentException
     *             if the displays are not the same.
     */
    public DisCoGraphRenderingCanvas(final Composite parent, final DisCoGraphRenderingConfigurator configurator) {
        super(parent, SWT.NONE);

        if (!configurator.getDisplay().equals(parent.getDisplay())) {
            throw new IllegalArgumentException("configurator and parent composite must use the same display.");
        }

        addPaintListener(this);
        graphRenderer = new DisCoGraphRenderer(configurator);

        // SUPPRESS CHECKSTYLE NEXT MagicNumber
        backgroundColor = new Color(parent.getDisplay(), 255, 255, 255);
        setBackground(backgroundColor);
    }

    @Override
    public void dispose() {
        graphRenderer.dispose();
        backgroundColor.dispose();
    }

    /**
     * Sets the given layout graph as the painted graph.
     * 
     * @param thelayoutGraph
     *            layout graph to be painted
     */
    public void setLayoutGraph(final ElkNode thelayoutGraph) {
        // set new size values for the canvas
        if (thelayoutGraph != null) {
            //KShapeLayout shapeLayout = thelayoutGraph.getData(KShapeLayout.class);
            setSize(new Point((int) ((thelayoutGraph.getWidth() + 1) * graphRenderer.getScale()),
                    (int) ((thelayoutGraph.getHeight() + 1) * graphRenderer.getScale())));
        }

        this.layoutGraph = thelayoutGraph;
        redraw();
    }

    /**
     * Returns the currently painted layout graph.
     * 
     * @return the painted layout graph
     */
    public ElkNode getLayoutGraph() {
        return layoutGraph;
    }

    // changed by mic
    /**
     * Returns the renderer used for painting.
     * 
     * @return the renderer
     */
    public DisCoGraphRenderer getRenderer() {
        return graphRenderer;
    }

    /**
     * This method is called when the canvas is requested to paint.
     * 
     * @param event
     *            paint event
     */
    @Override
    public void paintControl(final PaintEvent event) {
        if (layoutGraph != null) {
            Rectangle area = new Rectangle(event.x, event.y, event.width, event.height);
            // reset paint information
            graphRenderer.markDirty(area);
            // paint the top layout nodes with their children
            graphRenderer.render(layoutGraph, event.gc, area);

        }
    }

}
