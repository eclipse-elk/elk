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

import java.awt.geom.Path2D;
import java.awt.geom.Rectangle2D;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.eclipse.elk.alg.common.polyomino.structures.Direction;
import org.eclipse.elk.alg.common.polyomino.structures.Polyomino;
import org.eclipse.elk.alg.common.utils.UniqueTriple;
import org.eclipse.elk.alg.disco.graph.DCComponent;
import org.eclipse.elk.alg.disco.graph.DCDirection;
import org.eclipse.elk.alg.disco.graph.DCElement;
import org.eclipse.elk.alg.disco.graph.DCExtension;
import org.eclipse.elk.alg.disco.graph.DCGraph;
import org.eclipse.elk.alg.disco.options.DisCoOptions;
import org.eclipse.elk.alg.disco.structures.DCPolyomino;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.options.EdgeType;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.Triple;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.ElkShape;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.swt.SWT;
import org.eclipse.swt.graphics.GC;
import org.eclipse.swt.graphics.Rectangle;

import com.google.common.collect.Sets;

/**
 * Utility class that is able to render an ELK graph instance plus its
 * corresponding DCGraph and polyominoes, if these exist. This is primarily a
 * debug tool.
 * 
 * This class started as a copy of GraphRenderer.java, commit 71bb8c2f542,
 * 2016-05-07. Changes are made by mic and commented accordingly.
 */
public class DisCoGraphRenderer {

    /** default length of edge arrows. */
    private static final double ARROW_LENGTH = 8.0f;
    /** default width of edge arrows. */
    private static final double ARROW_WIDTH = 7.0f;
    /** the minimal font height for displaying labels. */
    private static final int MIN_FONT_HEIGHT = 3;

    /** mapping of each layout graph element to its computed bounds. */
    private final Map<Object, PaintRectangle> boundsMap = new LinkedHashMap<Object, PaintRectangle>();
    /** by mic: mapping each DCElement to its polygonal representation. */
    private final Map<Object, PaintPolygon> polygonMap = new LinkedHashMap<Object, PaintPolygon>();
    private final Map<Object, PaintRectangle> dcExtensionMap = new LinkedHashMap<Object, PaintRectangle>();
    /** by mic: mapping each polyomino cell to its computed bounds. */
    private final Map<Object, PaintRectangle> polyominoMap = new LinkedHashMap<Object, PaintRectangle>();
    private final Map<Object, PaintRectangle> polyominoCenterMap = new LinkedHashMap<Object, PaintRectangle>();
    private final Map<Object, PaintRectangle> polyominoExtensionMap = new LinkedHashMap<Object, PaintRectangle>();

    private ConfigState state;
    private DisCoGraphRenderingFillPatterns patterns;

    /**
     * configurator used to configure what the drawing looks like. Replaced by mic.
     */
    private DisCoGraphRenderingConfigurator configurator;
    /** the scale factor for all coordinates. */
    private double scale = 1.0;
    /** the base offset for all coordinates. */
    private KVector baseOffset = new KVector();
    /**
     * The most recently drawn graph; if the next graph to be drawn is different, we
     * will automatically flush our cache.
     */
    private ElkNode mostRecentlyDrawnGraph = null;

    /**
     * by mic: switches for turning certain representations of the graph on/off.
     */
    private boolean paintElkGraph = true;
    private boolean paintDCGraph = true;
    private boolean paintPolys = true;

    private Boolean savePaintDCGraph;
    private Boolean savePaintPolys;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialization / Cleanup

    /**
     * Creates new renderer that uses the given configurator to define the way
     * drawings look and a default scaling of 1.0.
     * 
     * @param configurator
     *            the rendering configurator.
     */
    public DisCoGraphRenderer(final DisCoGraphRenderingConfigurator configurator) {
        this(configurator, new ConfigState().getScale());
    }

    /**
     * Creates new renderer that uses the given configurator to define the way
     * drawings look.
     * 
     * @param configurator
     *            the rendering configurator.
     * @param scale
     *            the scaling factor to draw graphs with.
     */
    public DisCoGraphRenderer(final DisCoGraphRenderingConfigurator configurator, final double scale) {
        this.scale = scale;
        this.configurator = configurator;
        this.configurator.initialize(scale);
        this.state = new ConfigState();

        this.patterns = new DisCoGraphRenderingFillPatterns(configurator, scale * state.getPatternScale());
    }

    /**
     * Dispose all created resources such as colors and fonts.
     */
    public void dispose() {
        flushCache();
        mostRecentlyDrawnGraph = null;
        configurator.dispose();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Getters / Setters

    /**
     * Returns the base offset to be added to all element coordinates when drawing
     * them.
     */
    public KVector getBaseOffset() {
        return baseOffset;
    }

    /**
     * graphics.setClipping(rect.x, rect.y, rect.width, rect.height); Sets the base
     * offset to be added to all element coordinates when drawing them.
     * 
     * @param baseOffset
     *            the new base offset.
     */
    public void setBaseOffset(final KVector baseOffset) {
        if (baseOffset == null) {
            this.baseOffset = new KVector();
        } else {
            this.baseOffset = baseOffset;
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Drawing Control

    /**
     * Mark all objects in the given area as dirty.
     * 
     * @param area
     *            the area to mark as dirty, or {@code null} if all objects shall be
     *            marked
     */
    public void markDirty(final Rectangle area) {
        for (PaintRectangle rectangle : boundsMap.values()) {
            if (area == null || rectangle.intersects(area)) {
                rectangle.painted = false;
            }
        }
        for (PaintPolygon polygon : polygonMap.values()) {
            if (area == null || polygon.intersects(area)) {
                polygon.painted = false;
            }
        }
        for (PaintRectangle rectangle : polyominoMap.values()) {
            if (area == null || rectangle.intersects(area)) {
                rectangle.painted = false;
            }
        }
        for (PaintRectangle rectangle : dcExtensionMap.values()) {
            if (area == null || rectangle.intersects(area)) {
                rectangle.painted = false;
            }
        }
        for (PaintRectangle rectangle : polyominoExtensionMap.values()) {
            if (area == null || rectangle.intersects(area)) {
                rectangle.painted = false;
            }
        }
    }

    /**
     * Clear all internally cached data on painted graphs. Call this method if the
     * graph to be drawn next changes.
     */
    private void flushCache() {
        boundsMap.clear();
        polygonMap.clear();
        polyominoMap.clear();
        dcExtensionMap.clear();
        polyominoExtensionMap.clear();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Rendering Code

    /**
     * Paints the contained layout graph onto the given graphics object.
     * 
     * @param parentNode
     *            the parent node that shall be rendered
     * @param graphics
     *            the graphics context used to paint
     * @param area
     *            the area to fill
     */
    public void render(final ElkNode parentNode, final GC graphics, final Rectangle area) {
        int oldLineWidth = graphics.getLineWidth();
        graphics.setLineWidth((int) (scale * state.getThicknessScale()));

        /** by mic: the associated graph of connected components */
        DCGraph componentGraph = (DCGraph) parentNode.getProperty(DisCoOptions.DEBUG_DISCO_GRAPH);

        if (mostRecentlyDrawnGraph != parentNode) {
            flushCache();
            mostRecentlyDrawnGraph = parentNode;
        }

        // activate interpolation
        graphics.setInterpolation(SWT.HIGH);

        // determine an overall alpha value for nodes, depending on the
        // maximal node depth
        int maxDepth = Math.max(maxDepth(parentNode), 1);
        state.setDepth(maxDepth);
        // CHECKSTYLEOFF MagicNumber
        int nodeAlpha = 200 / maxDepth + 55;
        // by mic: all elements drawn will be quite transparent, so that the
        // user can see the different layers (graph, DCCGraph und polyomino
        // cells) stacked on top of each other.
        int fillingAlpha = 40 / maxDepth + 55;

        if (state.makeSolid()) {
            nodeAlpha = 0xff;
        }
        // CHECKSTYLEON MagicNumber

        if (componentGraph == null) {
            savePaintDCGraph = Boolean.valueOf(paintDCGraph);
            paintDCGraph = false;
            savePaintPolys = Boolean.valueOf(paintPolys);
            paintPolys = false;
        } else {
            if (savePaintDCGraph != null) {
                paintDCGraph = savePaintDCGraph.booleanValue();
                savePaintDCGraph = null;
            }
            if (savePaintPolys != null) {
                paintPolys = savePaintPolys.booleanValue();
                savePaintPolys = null;
            }
        }

        renderRecursively(parentNode, graphics, area, baseOffset, nodeAlpha, fillingAlpha, 0);

        if (componentGraph == null) {
            graphics.setForeground(configurator.getBlack());
            graphics.drawString("No layout algorithm for connected components has been used.", 0, 0, true);
        }

        graphics.setLineWidth(oldLineWidth);
    }

    private void renderRecursively(final ElkNode parent, final GC graphics, final Rectangle area, final KVector offset,
            final int nodeAlpha, final int fillingAlpha, final int levelNumber) {
        /** by mic: the associated graph of connected components */
        DCGraph componentGraph = (DCGraph) parent.getProperty(DisCoOptions.DEBUG_DISCO_GRAPH);

        int lowerLvl = state.getLowerLvl();
        int upperLvl = state.getUpperLvl();

        if ((paintElkGraph && levelNumber >= lowerLvl && levelNumber <= upperLvl)
                || (state.drawGhostParent() && levelNumber == 0)) {
            Set<ElkEdge> edgeSet = new HashSet<ElkEdge>();
            // render the nodes and ports
            if (lowerLvl > 0 && state.drawGhostParent() && levelNumber == 0) {
                renderNodeChildren(parent, graphics, area, offset, edgeSet, nodeAlpha, 0);
            } else {
                renderNodeChildren(parent, graphics, area, offset, edgeSet, nodeAlpha, fillingAlpha);
            }

            // render the edges
            for (ElkEdge edge : edgeSet) {
                renderEdge(parent, edge, graphics, area, offset, nodeAlpha);
            }
        }

        if (paintDCGraph && componentGraph != null && levelNumber >= lowerLvl && levelNumber <= upperLvl) {
            renderComponentGraph(parent, componentGraph, graphics, area, offset, nodeAlpha, fillingAlpha, levelNumber);
        }

        if (paintPolys && levelNumber >= lowerLvl && levelNumber <= upperLvl) {
            // by mic: low resolution polyominos generated from the graph of
            // connected components
            @SuppressWarnings("unchecked")
            List<DCPolyomino> polys = (List<DCPolyomino>) parent.getProperty(DisCoOptions.DEBUG_DISCO_POLYS);
            if (polys != null) {
                renderPolyominoes(parent, polys, graphics, area, offset, nodeAlpha, fillingAlpha, levelNumber);
            }
        }

        for (ElkNode child : parent.getChildren()) {
            // compute the offset required to make the children's
            KVector contentOffset = new KVector(child.getX() * getScale(), child.getY() * getScale());
            contentOffset.add(offset);

            renderRecursively(child, graphics, area, contentOffset, nodeAlpha, fillingAlpha, levelNumber + 1);
        }
    }

    // added by mic
    /**
     * Paints all elements of the DCGraph that fall into the given dirty area.
     * 
     * @param componentGraph
     *            The DCGraph to paint
     * @param graphics
     *            the graphics context used to paint
     * @param area
     *            dirty area that needs painting
     * @param offset
     *            offset to be added to relative coordinates
     * @param nodeAlpha
     *            alpha value for nodes
     * @param fillingAlpha
     *            alpha value for node fillings
     */
    private void renderComponentGraph(final ElkNode parent, final DCGraph componentGraph, final GC graphics,
            final Rectangle area, final KVector offset, final int nodeAlpha, final int fillingAlpha,
            final int levelNumber) {
        KVector boundaries = new KVector(parent.getWidth(), parent.getHeight());

        for (DCComponent comp : componentGraph.getComponents()) {

            for (DCElement el : comp.getElements()) {
                PaintPolygon poly = polygonMap.get(el);
                if (poly == null) {
                    poly = new PaintPolygon(el, offset, getScale());
                    polygonMap.put(el, poly);
                }

                if (!poly.painted && poly.intersects(area)) {
                    // paint this node
                    graphics.setAlpha(fillingAlpha);

                    if (configurator.getDCElementFillColor() != null) {
                        graphics.setBackground(configurator.getDCElementFillColor());
                        graphics.fillPolygon(poly.coordsScaledAndRounded);
                    }

                    graphics.setAlpha(nodeAlpha);

                    if (configurator.getDCElementBorderTextColor() != null) {
                        graphics.setForeground(configurator.getDCElementBorderTextColor());
                        graphics.drawPolygon(poly.coordsScaledAndRounded);
                        if (configurator.getNodeLabelFont() != null) {
                            graphics.setFont(configurator.getNodeLabelFont());
                        }
                        if (state.drawLabels()) {
                            String levelprefix = levelNumber + "_";
                            if (state.removeLvl()) {
                                levelprefix = "";
                            }
                            graphics.drawString(levelprefix + Integer.toString(comp.getId()),
                                    poly.coordsScaledAndRounded[0], poly.coordsScaledAndRounded[1], true);
                        }
                    }

                    poly.painted = true;

                }

                ElkRectangle bounds = el.getBounds();
                Rectangle2D.Double elementBounding = new Rectangle2D.Double(bounds.x, bounds.y, bounds.width,
                        bounds.height);
                double offsetX = elementBounding.getX() + el.getOffset().x;
                double offsetY = elementBounding.getY() + el.getOffset().y;
                elementBounding = new Rectangle2D.Double(offsetX, offsetY, elementBounding.getWidth(),
                        elementBounding.getHeight());

                for (DCExtension ext : el.getExtensions()) {

                    PaintRectangle rect = dcExtensionMap.get(ext);
                    if (rect == null) {
                        rect = new PaintRectangle(ext, elementBounding, boundaries, offset, getScale());
                        dcExtensionMap.put(ext, rect);
                    }

                    if (!rect.painted && rect.intersects(area)) {
                        graphics.setAlpha(fillingAlpha);

                        // CHECKSTYLEOFF MagicNumber
                        if (configurator.getDCElementExternalFillColor() != null) {
                            graphics.setBackgroundPattern(
                                    patterns.getDCExtensionPattern(state.makeSolid() ? 255 : fillingAlpha));
                            graphics.fillRectangle(rect.x, rect.y, rect.width, rect.height);
                        }
                        // CHECKSTYLEON MagicNumber

                        graphics.setAlpha(nodeAlpha);

                        if (configurator.getNodeBorderColor() != null) {
                            graphics.setForeground(configurator.getDCElementExternalBorderTextColor());
                            graphics.drawRectangle(rect.x, rect.y, rect.width, rect.height);
                        }

                        if (configurator.getDCElementExternalBorderTextColor() != null) {
                            graphics.setForeground(configurator.getDCElementExternalBorderTextColor());
                            if (configurator.getNodeLabelFont() != null) {
                                graphics.setFont(configurator.getNodeLabelFont());
                            }
                            if (state.drawLabels()) {
                                String levelprefix = levelNumber + "_";
                                if (state.removeLvl()) {
                                    levelprefix = "";
                                }
                                graphics.drawString(levelprefix + Integer.toString(comp.getId()), rect.x, rect.y, true);
                            }
                        }

                        rect.painted = true;

                    }

                }
            }
        }
    }

    // added by mic
    /**
     * Paints all polyominoes that fall into the given dirty area.
     * 
     * @param polys
     *            the polyominoes to paint
     * @param graphics
     *            the graphics context used to paint
     * @param area
     *            dirty area that needs painting
     * @param offset
     *            offset to be added to relative coordinates
     * @param nodeAlpha
     *            alpha value for nodes
     * @param fillingAlpha
     *            alpha value for node fillings
     */
    private void renderPolyominoes(final ElkNode parent, final List<DCPolyomino> polys, final GC graphics,
            final Rectangle area, final KVector offset, final int nodeAlpha, final int fillingAlpha,
            final int levelNumber) {

        for (DCPolyomino poly : polys) {
            KVector polyCorner = poly.getMinCornerOnCanvas();
            double topLeftCornerX = polyCorner.x;
            double topLeftCornerY = polyCorner.y;
            double cellSizeX = poly.getCellSizeX();
            double cellSizeY = poly.getCellSizeY();
            double bottomRightCornerX = topLeftCornerX + cellSizeX * poly.getWidth();
            double bottomRightCornerY = topLeftCornerY + cellSizeY * poly.getHeight();

            for (int x = 0; x < poly.getWidth(); x++) {
                for (int y = 0; y < poly.getHeight(); y++) {
                    Triple<Polyomino, Integer, Integer> polyoCell = new Triple<Polyomino, Integer, Integer>(poly, x, y);
                    PaintRectangle rect = polyominoMap.get(polyoCell);
                    if (rect == null) {
                        rect = new PaintRectangle(topLeftCornerX, topLeftCornerY, cellSizeX, cellSizeY, x, y, offset,
                                getScale());
                        polyominoMap.put(polyoCell, rect);
                    }

                    if (!rect.painted && rect.intersects(area)) {
                        // paint this node
                        graphics.setAlpha(fillingAlpha);

                        if (configurator.getPolyominoFillColor() != null && poly.isBlocked(x, y)) {
                            graphics.setBackground(configurator.getPolyominoFillColor());
                            graphics.fillRectangle(rect.x, rect.y, rect.width, rect.height);
                            if (configurator.getNodeLabelFont() != null) {
                                graphics.setFont(configurator.getNodeLabelFont());
                            }
                            graphics.setAlpha(nodeAlpha);
                            graphics.setForeground(configurator.getPolyominoBorderTextColor());
                            if (state.drawLabels()) {
                                String levelprefix = levelNumber + "_";
                                if (state.removeLvl()) {
                                    levelprefix = "";
                                }
                                graphics.drawString(levelprefix + Integer.toString(poly.getId()), rect.x, rect.y, true);
                            }
                        }

                        if (configurator.getPolyominoFillColor() != null && poly.isWeaklyBlocked(x, y)) {
                            // CHECKSTYLEOFF MagicNumber
                            graphics.setBackgroundPattern(
                                    patterns.getPolyominoExtensionPattern(state.makeSolid() ? 255 : fillingAlpha));
                            // CHECKSTYLEON MagicNumber
                            graphics.fillRectangle(rect.x, rect.y, rect.width, rect.height);
                            graphics.setAlpha(nodeAlpha);
                            graphics.setForeground(configurator.getPolyominoWeaklyBlockedBorderTextColor());
                            graphics.drawRectangle(rect.x, rect.y, rect.width, rect.height);

                            if (configurator.getNodeLabelFont() != null) {
                                graphics.setFont(configurator.getNodeLabelFont());
                            }
                            graphics.setAlpha(nodeAlpha);
                            if (state.drawLabels()) {
                                String levelprefix = levelNumber + "_";
                                if (state.removeLvl()) {
                                    levelprefix = "";
                                }
                                graphics.drawString(levelprefix + Integer.toString(poly.getId()), rect.x, rect.y, true);
                            }

                        }

                        graphics.setAlpha(nodeAlpha);

                        if (configurator.getPolyominoBorderTextColor() != null && !poly.isWeaklyBlocked(x, y)) {
                            graphics.setForeground(configurator.getPolyominoBorderTextColor());

                            if (state.drawPolyLinesBlack()) {
                                graphics.setForeground(configurator.getBlack());
                            }

                            graphics.drawRectangle(rect.x, rect.y, rect.width, rect.height);
                        }

                        rect.painted = true;
                    }
                }
            }

            if (state.markTheCenter()) {
                int x = poly.getCenterX();
                int y = poly.getCenterY();

                PaintRectangle centerRect = polyominoCenterMap.get(poly);
                if (centerRect == null) {
                    centerRect = new PaintRectangle(topLeftCornerX, topLeftCornerY, cellSizeX, cellSizeY, x, y, offset,
                            getScale());
                    polyominoCenterMap.put(poly, centerRect);
                }

                if (configurator.getPolyominoFillColor() != null) {
                    // CHECKSTYLEOFF MagicNumber
                    graphics.setBackgroundPattern(
                            patterns.getPolyominoCenterPattern(state.makeSolid() ? 255 : fillingAlpha));
                    // CHECKSTYLEON MagicNumber
                    graphics.fillRectangle(centerRect.x, centerRect.y, centerRect.width, centerRect.height);
                    graphics.setAlpha(nodeAlpha);

                }
            }

            for (UniqueTriple<Direction, Integer, Integer> ext : poly.getPolyominoExtensions()) {
                PaintRectangle rect = polyominoExtensionMap.get(ext);
                if (rect == null) {

                    Direction dir = ext.getFirst();
                    int extStart = ext.getSecond();
                    int extEnd = ext.getThird();

                    switch (dir) {
                    case NORTH:
                        rect = new PaintRectangle(topLeftCornerX + extStart * cellSizeX, 0,
                                (extEnd - extStart + 1) * cellSizeX, topLeftCornerY, offset, getScale());
                        break;

                    case SOUTH:
                        rect = new PaintRectangle(topLeftCornerX + extStart * cellSizeX, bottomRightCornerY,
                                (extEnd - extStart + 1) * cellSizeX, parent.getHeight() - bottomRightCornerY, offset,
                                getScale());
                        break;
                    case WEST:
                        rect = new PaintRectangle(0, topLeftCornerY + extStart * cellSizeY, topLeftCornerX,
                                (extEnd - extStart + 1) * cellSizeY, offset, getScale());
                        break;

                    default: // EAST
                        rect = new PaintRectangle(bottomRightCornerX, topLeftCornerY + extStart * cellSizeY,
                                parent.getWidth() - bottomRightCornerX, (extEnd - extStart + 1) * cellSizeY, offset,
                                getScale());
                        break;

                    }
                    // Maybe not unique in extreme cases where extensions
                    // are very close to each other, but it doesn't matter
                    // for the drawing
                    polyominoExtensionMap.put(ext, rect);

                }

                graphics.setAlpha(fillingAlpha);
                // CHECKSTYLEOFF MagicNumber
                graphics.setBackgroundPattern(
                        patterns.getPolyominoExtensionPattern(state.makeSolid() ? 255 : fillingAlpha));
                // CHECKSTYLEON MagicNumber
                graphics.fillRectangle(rect.x, rect.y, rect.width, rect.height);
                graphics.setAlpha(nodeAlpha);
                graphics.setForeground(configurator.getPolyominoWeaklyBlockedBorderTextColor());
                graphics.drawRectangle(rect.x, rect.y, rect.width, rect.height);

                if (configurator.getNodeLabelFont() != null) {
                    graphics.setFont(configurator.getNodeLabelFont());
                }
                graphics.setAlpha(nodeAlpha);
                if (state.drawLabels()) {
                    String levelprefix = levelNumber + "_";
                    if (state.removeLvl()) {
                        levelprefix = "";
                    }
                    graphics.drawString(levelprefix + Integer.toString(poly.getId()), rect.x, rect.y, true);
                }
            }
        }
    }

    // filling alpha added by mic
    /**
     * Paints all children of the given parent node that fall into the given dirty
     * area.
     * 
     * @param parent
     *            the node whose children to paint
     * @param graphics
     *            the graphics context used to paint
     * @param area
     *            dirty area that needs painting
     * @param offset
     *            offset to be added to relative coordinates
     * @param edgeSet
     *            set to be filled with edges that are found on the way
     * @param nodeAlpha
     *            alpha value for nodes
     * @param fillingAlpha
     *            alpha value for node fillings
     */
    private void renderNodeChildren(final ElkNode parent, final GC graphics, final Rectangle area, final KVector offset,
            final Set<ElkEdge> edgeSet, final int nodeAlpha, final int fillingAlpha) {

        for (ElkNode child : parent.getChildren()) {
            PaintRectangle rect = boundsMap.get(child);
            if (rect == null) {
                rect = new PaintRectangle(child, offset, getScale());
                boundsMap.put(child, rect);
            }
            KVector childOffset = new KVector(rect.x, rect.y);

            // render the child node and its content
            if (!rect.painted && rect.intersects(area)) {
                // paint this node
                graphics.setAlpha(fillingAlpha);

                if (configurator.getNodeFillColor() != null) {
                    graphics.setBackground(configurator.getNodeFillColor());
                    graphics.fillRectangle(rect.x, rect.y, rect.width, rect.height);
                }

                graphics.setAlpha(nodeAlpha);

                if (configurator.getNodeBorderColor() != null) {
                    graphics.setForeground(configurator.getNodeBorderColor());
                    graphics.drawRectangle(rect.x, rect.y, rect.width, rect.height);
                }

                rect.painted = true;
            }

            // render node labels
            if (configurator.getNodeLabelFont() != null) {
                graphics.setFont(configurator.getNodeLabelFont());
                for (ElkLabel label : child.getLabels()) {
                    renderLabel(label, graphics, area, childOffset, nodeAlpha);
                }
            }

            // render ports
            for (ElkPort port : child.getPorts()) {
                renderPort(port, graphics, area, childOffset, nodeAlpha);
            }

            // store all incident edges to render them later
            edgeSet.addAll(Sets.newHashSet(ElkGraphUtil.allIncidentEdges(child)));
        }
    }

    /**
     * Paints a label for the given dirty area. Expects the correct font to be
     * already set.
     * 
     * @param label
     *            the label to paint
     * @param graphics
     *            the graphics context used to paint
     * @param area
     *            dirty area that needs painting
     * @param offset
     *            offset to be added to relative coordinates
     * @param labelAlpha
     *            alpha value for labels
     */
    private void renderLabel(final ElkLabel label, final GC graphics, final Rectangle area, final KVector offset,
            final int labelAlpha) {

        if (graphics.getFont().getFontData()[0].getHeight() >= MIN_FONT_HEIGHT) {
            PaintRectangle rect = boundsMap.get(label);
            if (rect == null) {
                rect = new PaintRectangle(label, offset, getScale());
                boundsMap.put(label, rect);
            }

            if (!rect.painted && rect.intersects(area)) {
                // render the border and filling
                graphics.setAlpha(labelAlpha);

                if (configurator.getLabelFillColor() != null) {
                    graphics.setBackground(configurator.getLabelFillColor());
                    graphics.fillRectangle(rect.x, rect.y, rect.width, rect.height);
                }

                if (configurator.getLabelBorderColor() != null) {
                    graphics.setForeground(configurator.getLabelBorderColor());
                    graphics.drawRectangle(rect.x, rect.y, rect.width, rect.height);
                }

                // render the text
                String text = label.getText();
                if (text != null && text.length() > 0) {
                    graphics.setForeground(configurator.getLabelTextColor());

                    Rectangle oldClip = graphics.getClipping();

                    if (state.drawLabels()) {
                        graphics.drawString(text, rect.x, rect.y, true);
                    }
                    graphics.setClipping(oldClip);
                }

                rect.painted = true;
            }
        }
    }

    /**
     * Paints a port for the given dirty area.
     * 
     * @param port
     *            the port to paint
     * @param graphics
     *            the graphics context used to paint
     * @param area
     *            dirty area that needs painting
     * @param offset
     *            offset to be added to relative coordinates
     * @param labelAlpha
     *            alpha value for labels
     */
    private void renderPort(final ElkPort port, final GC graphics, final Rectangle area, final KVector offset,
            final int labelAlpha) {

        PaintRectangle rect = boundsMap.get(port);
        if (rect == null) {
            rect = new PaintRectangle(port, offset, getScale());
            boundsMap.put(port, rect);
        }

        if (!rect.painted && rect.intersects(area)) {
            // CHECKSTYLEOFF MagicNumber
            graphics.setAlpha(255);
            // CHECKSTYLEON MagicNumber

            if (configurator.getPortFillColor() != null) {
                graphics.setBackground(configurator.getPortFillColor());
                graphics.fillRectangle(rect.x, rect.y, rect.width, rect.height);
            }

            if (configurator.getPortBorderColor() != null) {
                graphics.setForeground(configurator.getPortBorderColor());
                graphics.drawRectangle(rect.x, rect.y, rect.width, rect.height);
            }

            rect.painted = true;
        }

        // paint port labels
        if (configurator.getPortLabelFont() != null) {
            graphics.setFont(configurator.getPortLabelFont());
            KVector portOffset = new KVector(rect.x, rect.y);
            for (ElkLabel label : port.getLabels()) {
                renderLabel(label, graphics, area, portOffset, labelAlpha);
            }
        }
    }

    /**
     * Paints an edge for the given dirty area.
     * 
     * @param graph
     *            the top-level node of the graph
     * @param edge
     *            the edge to paint
     * @param graphics
     *            the graphics context used to paint
     * @param area
     *            dirty area that needs painting
     * @param labelAlpha
     *            alpha value for labels
     */
    private void renderEdge(final ElkNode graph, final ElkEdge edge, final GC graphics, final Rectangle area,
            final KVector edgeBaseOffset, final int labelAlpha) {

        if (configurator.getEdgeColor() == null) {
            return;
        }

        if (!ElkGraphUtil.isDescendant(ElkGraphUtil.getSourceNode(edge), graph)
                || !ElkGraphUtil.isDescendant(ElkGraphUtil.getTargetNode(edge), graph)) {

            // the edge points to some node outside of the rendered subgraph
            return;
        }

        // calculate an offset for edge coordinates
        ElkNode parent = ElkGraphUtil.getSourceNode(edge);
        if (!ElkGraphUtil.isDescendant(ElkGraphUtil.getTargetNode(edge), parent)) {
            parent = parent.getParent();
        }
        ElkNode node = parent;
        KVector offset = new KVector().add(edgeBaseOffset);
        while (node != graph) {
            offset.add(node.getX() * getScale(), node.getY() * getScale());
            node = node.getParent();
        }

        // ElkEdgeLayout edgeLayout = edge.getData(ElkEdgeLayout.class);
        PaintRectangle rect = boundsMap.get(edge);

        if (rect == null) {
            rect = new PaintRectangle(edge, offset, getScale());
            boundsMap.put(edge, rect);
        }

        if (!rect.painted && rect.intersects(area)) {
            // CHECKSTYLEOFF MagicNumber
            graphics.setAlpha(255);
            // CHECKSTYLEON MagicNumber

            // The background color is required to fill the arrow of directed
            // edges
            graphics.setForeground(configurator.getEdgeColor());
            graphics.setBackground(configurator.getEdgeColor());

            ElkEdgeSection edgeSection = ElkGraphUtil.firstEdgeSection(edge, false, false);
            KVectorChain bendPoints = ElkUtil.createVectorChain(edgeSection);
            if (edge.getProperty(CoreOptions.EDGE_ROUTING) == EdgeRouting.SPLINES) {
                bendPoints = ElkMath.approximateBezierSpline(bendPoints);
            }
            bendPoints.scale(getScale()).offset(offset);
            KVector point1 = bendPoints.getFirst();
            for (KVector point2 : bendPoints) {
                graphics.drawLine((int) Math.round(point1.x), (int) Math.round(point1.y), (int) Math.round(point2.x),
                        (int) Math.round(point2.y));
                point1 = point2;
            }

            if (edge.getProperty(CoreOptions.EDGE_TYPE) != EdgeType.UNDIRECTED) {
                // draw an arrow at the last segment of the connection
                int[] arrowPoly = makeArrow(bendPoints.get(bendPoints.size() - 2), bendPoints.getLast());
                if (arrowPoly != null) {
                    graphics.fillPolygon(arrowPoly);
                }
            }

            rect.painted = true;
        }

        // paint junction points
        KVectorChain vc = edge.getProperty(CoreOptions.JUNCTION_POINTS);
        if (vc != null) {
            for (KVector v : vc) {
                KVector center = v.clone().scale(getScale()).add(offset).sub(2, 2);
                // CHECKSTYLEOFF MagicNumber
                graphics.fillOval((int) center.x, (int) center.y, 6, 6);
                // CHECKSTYLEON MagicNumber
            }
        }

        // paint the edge labels
        if (configurator.getEdgeLabelFont() != null) {
            graphics.setFont(configurator.getEdgeLabelFont());
            for (ElkLabel label : edge.getLabels()) {
                renderLabel(label, graphics, area, offset, labelAlpha);
            }
        }
    }

    /**
     * Constructs a polygon that forms an arrow.
     * 
     * @param point1
     *            source point
     * @param point2
     *            target point
     * @return array of coordinates for the arrow polygon, or null if the given
     *         source and target points are equal
     */
    private int[] makeArrow(final KVector point1, final KVector point2) {
        if (!(point1.x == point2.x && point1.y == point2.y) && ARROW_WIDTH * getScale() >= 2) {
            // CHECKSTYLEOFF MagicNumber
            int[] arrow = new int[6];
            arrow[0] = (int) Math.round(point2.x);
            arrow[1] = (int) Math.round(point2.y);

            double vectX = point1.x - point2.x;
            double vectY = point1.y - point2.y;
            double length = Math.sqrt(vectX * vectX + vectY * vectY);
            double normX = vectX / length;
            double normY = vectY / length;
            double neckX = point2.x + ARROW_LENGTH * normX * getScale();
            double neckY = point2.y + ARROW_LENGTH * normY * getScale();
            double orthX = normY * ARROW_WIDTH / 2 * getScale();
            double orthY = -normX * ARROW_WIDTH / 2 * getScale();

            arrow[2] = (int) Math.round(neckX + orthX);
            arrow[3] = (int) Math.round(neckY + orthY);
            arrow[4] = (int) Math.round(neckX - orthX);
            arrow[5] = (int) Math.round(neckY - orthY);
            // CHECKSTYLEON MagicNumber
            return arrow;
        } else {
            return null;
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Utility Things

    /**
     * Determine the maximal depth of the given graph.
     * 
     * @param parent
     *            the parent node of the graph
     * @return the maximal depth of contained nodes
     */
    private int maxDepth(final ElkNode parent) {
        int maxDepth = 0;
        for (ElkNode child : parent.getChildren()) {
            int depth = maxDepth(child) + 1;
            if (depth > maxDepth) {
                maxDepth = depth;
            }
        }
        return maxDepth;
    }

    /**
     * Method for switching aspects of the drawing on and off. Commands:</br>
     * 
     * <ul>
     * <li>"kgraph": toggles the visibility of the original layout graph.</li>
     * <li>"dcgraph": toggles the visibility of the associated DCGraph.</li>
     * <li>"polyominoes": toggles the visibility of the polyominoes.</li>
     * </ul>
     * 
     * @param command
     *            Determines the aspect of the drawing whose visibility is to be
     *            toggled.
     * @param visible
     *            true: show drawing, false: hide drawing
     */
    public void setVisible(final String command, final boolean visible) {
        switch (command) {
        case "elkgraph":
            this.paintElkGraph = visible;
            break;
        case "dcgraph":
            this.paintDCGraph = visible;
            break;
        case "polyominoes":
            this.paintPolys = visible;
        default:
            break;
        }
    }

    /**
     * Return the configuration object.
     * 
     * @return the configuration state
     */
    ConfigState getState() {
        return state;
    }

    /**
     * Return the scaling factor for all coordinates.
     * 
     * @return the scaling factor
     */
    public double getScale() {
        return scale;
    }

    /**
     * Rectangle class used to mark painted objects.
     */
    private static class PaintRectangle {
        private int x, y, width, height;
        private boolean painted = false;

        PaintRectangle(final double x, final double y, final double width, final double height, final KVector offset,
                final double scale) {
            this.x = (int) Math.round(x * scale + offset.x);
            this.y = (int) Math.round(y * scale + offset.y);
            this.width = Math.max((int) Math.round(width * scale), 1);
            this.height = Math.max((int) Math.round(height * scale), 1);
        }

        /**
         * Creates a paint rectangle from a shape layout object.
         * 
         * @param shapeLayout
         *            shape layout from which values shall be taken
         * @param offset
         *            offset to be added to coordinate values
         */
        PaintRectangle(final ElkShape shapeLayout, final KVector offset, final double scale) {
            this.x = (int) Math.round(shapeLayout.getX() * scale + offset.x);
            this.y = (int) Math.round(shapeLayout.getY() * scale + offset.y);
            this.width = Math.max((int) Math.round(shapeLayout.getWidth() * scale), 1);
            this.height = Math.max((int) Math.round(shapeLayout.getHeight() * scale), 1);
        }

        /**
         * Creates a paint rectangle from an edge layout object.
         * 
         * @param edgeLayout
         *            edge layout from which the values shall be determined
         * @param offset
         *            offset to be added to coordinate values
         * @param scale
         *            the scale to apply to all coordinates
         */
        PaintRectangle(final ElkEdge edge, final KVector offset, final double scale) {
            ElkNode source = ElkGraphUtil.getSourceNode(edge);
            ElkNode target = ElkGraphUtil.getTargetNode(edge);
            ElkEdgeSection edgeSection = ElkGraphUtil.firstEdgeSection(edge, false, false);
            KVectorChain bendPoints = ElkUtil.createVectorChain(edgeSection);
            double minX = source.getX(), minY = source.getY();
            double maxX = minX, maxY = minY;
            for (KVector point : bendPoints) {
                minX = Math.min(minX, point.x);
                minY = Math.min(minY, point.y);
                maxX = Math.max(maxX, point.x);
                maxY = Math.max(maxY, point.y);
            }
            minX = Math.min(minX, target.getX());
            minY = Math.min(minY, target.getY());
            maxX = Math.max(maxX, target.getX());
            maxY = Math.max(maxY, target.getY());
            this.x = (int) Math.round(minX * scale + offset.x);
            this.y = (int) Math.round(minY * scale + offset.y);
            this.width = (int) Math.round((maxX - minX) * scale);
            this.height = (int) Math.round((maxY - minY) * scale);
        }

        /**
         * Creates a paint rectangle from a polyomino cell.
         * 
         * @param polyCornerX
         *            Upper left corner of the polyomino
         * @param polyCornerY
         *            Upper left corner of the polyomino
         * @param cellSizeX
         *            Width of cell
         * @param cellSizeY
         *            Height of cell
         * @param cellX
         *            position of the single polyomino cell in the polyomino's grid
         * @param cellY
         *            position of the single polyomino cell in the polyomino's grid
         * @param offset
         *            offset to be added to coordinate values
         * @param scale
         *            the scale to apply to all coordinates
         */
        PaintRectangle(final double polyCornerX, final double polyCornerY, final double cellSizeX,
                final double cellSizeY, final int cellX, final int cellY, final KVector offset, final double scale) {
            this.x = (int) Math.round((polyCornerX + cellX * cellSizeX) * scale + offset.x);
            this.y = (int) Math.round((polyCornerY + cellY * cellSizeY) * scale + offset.y);
            this.width = Math.max((int) Math.round(cellSizeX * scale), 1);
            this.height = Math.max((int) Math.round(cellSizeY * scale), 1);
        }

        /**
         * Creates a paint rectangle from an extension.
         * 
         * @param ext
         *            the extension
         * @param elementBounding
         *            bounding rectangle of the node the extension belongs to
         * @param boundaries
         *            width and height of surrounding parent node
         * @param offset
         *            offset to be added to coordinate values
         * @param scale
         *            the scale to apply to all coordinates
         */
        PaintRectangle(final DCExtension ext, final Rectangle2D elementBounding, final KVector boundaries,
                final KVector offset, final double scale) {
            DCDirection dir = ext.getDirection();
            KVector offsetFromDcElement = ext.getOffset();
            double widthExt = ext.getWidth();

            double xD, yD, widthD, heightD;

            switch (dir) {
            case NORTH:
                xD = elementBounding.getX() + offsetFromDcElement.x;
                yD = 0;
                widthD = widthExt;
                heightD = elementBounding.getY() + offsetFromDcElement.y;
                break;
            case SOUTH:
                xD = elementBounding.getX() + offsetFromDcElement.x;
                yD = elementBounding.getY() + offsetFromDcElement.y;
                widthD = widthExt;
                heightD = boundaries.y - yD;
                break;
            case WEST:
                xD = 0;
                yD = elementBounding.getY() + offsetFromDcElement.y;
                widthD = elementBounding.getX() + offsetFromDcElement.x - xD;
                heightD = widthExt;
                break;
            default: // EAST
                xD = elementBounding.getX() + offsetFromDcElement.x;
                yD = elementBounding.getY() + offsetFromDcElement.y;
                widthD = boundaries.x - xD;
                heightD = widthExt;
            }

            this.x = (int) Math.round(xD * scale + offset.x);
            this.y = (int) Math.round(yD * scale + offset.y);
            this.width = Math.max((int) Math.round(widthD * scale), 1);
            this.height = Math.max((int) Math.round(heightD * scale), 1);
        }

        /**
         * Determines whether the given rectangle intersects with the receiver.
         * 
         * @param other
         *            the rectangle to check for intersection
         * @return true if the other rectangle intersects with this one
         */
        public boolean intersects(final Rectangle other) {
            return (other.x < this.x + this.width) && (other.y < this.y + this.height)
                    && (other.x + other.width > this.x) && (other.y + other.height > this.y);
        }
    }

    /**
     * Polygon class to represent {@link DCElement DCElements} to be drawn and used
     * to mark painted shapes.
     * 
     * @author mic
     *
     */
    private static class PaintPolygon {
        /**
         * Coordinates of the DCElement translated to the coordinate system of the
         * canvas. Used to paint the shape.
         */
        private int[] coordsScaledAndRounded;
        /** Indicates whether the polygon has been painted yet. */
        private boolean painted = false;
        /**
         * Scaled coordinates as a {@link Path2D}, for using its provided method for
         * testing for intersections.
         */
        private Path2D.Float pathRepresentation = new Path2D.Float();

        /**
         * Creates a paint polygon from a {@link DCElement}.
         * 
         * @param elem
         *            {@link DCElement} from which the values shall be determined
         * @param offset
         *            offset to be added to coordinate values
         * @param scale
         *            the scale to apply to all coordinates
         */
        PaintPolygon(final DCElement elem, final KVector offset, final double scale) {
            double[] coords = elem.getCoords();
            KVector componentOffset = elem.getOffset();
            int len = coords.length;
            double offX = offset.x;
            double offY = offset.y;
            coordsScaledAndRounded = new int[len];
            for (int i = 0; i < len; i += 2) {
                coordsScaledAndRounded[i] = (int) Math.round((coords[i] + componentOffset.x) * scale + offX);
                coordsScaledAndRounded[i + 1] = (int) Math.round((coords[i + 1] + componentOffset.y) * scale + offY);
            }
            pathRepresentation.moveTo(coordsScaledAndRounded[0], coordsScaledAndRounded[1]);
            for (int j = 2; j < len; j += 2) {
                pathRepresentation.lineTo(coordsScaledAndRounded[j], coordsScaledAndRounded[j + 1]);
            }
            pathRepresentation.closePath();
        }

        /**
         * Determines whether the given rectangle intersects with the receiver.
         * 
         * @param other
         *            the rectangle to check for intersection
         * @return true if the other rectangle intersects with this polygon
         */
        public boolean intersects(final Rectangle other) {
            return pathRepresentation.intersects(new Rectangle2D.Float(other.x, other.y, other.width, other.height));
        }
    }
}
