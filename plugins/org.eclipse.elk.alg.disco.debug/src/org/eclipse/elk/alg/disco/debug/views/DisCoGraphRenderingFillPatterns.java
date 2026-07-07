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

import java.util.Arrays;

import org.eclipse.elk.alg.disco.graph.DCExtension;
import org.eclipse.swt.graphics.Color;
import org.eclipse.swt.graphics.GC;
import org.eclipse.swt.graphics.Image;
import org.eclipse.swt.graphics.ImageData;
import org.eclipse.swt.graphics.PaletteData;
import org.eclipse.swt.graphics.Pattern;
import org.eclipse.swt.widgets.Display;

/**
 * Class providing fill patterns for extensions.
 */
public class DisCoGraphRenderingFillPatterns {
    private DisCoGraphRenderingConfigurator config;
    private Display display;
    private int scale;

    /**
     * Creates a pattern creation object.
     * 
     * @param config
     *            configuration of the debug view
     * @param scale
     *            scaling factor for the pattern
     */
    public DisCoGraphRenderingFillPatterns(final DisCoGraphRenderingConfigurator config, final double scale) {
        this.config = config;
        this.display = config.getDisplay();
        this.scale = (int) scale;
    }

    /**
     * Get the pattern for a {@link DCExtension}.
     * 
     * @param alpha
     *            transparency
     * @return the pattern
     */
    public Pattern getDCExtensionPattern(final int alpha) {
        // This is hacky:
        // http://stackoverflow.com/questions/8270474/drawing-on-a-transparent-image-using-java-swt
        // but we want some totally transparent pixels
        // allocate an image data
        // CHECKSTYLEOFF MagicNumber
        ImageData imData = new ImageData(6 * scale, 6 * scale, 24, new PaletteData(0xff0000, 0x00ff00, 0x0000ff));
        imData.setAlpha(0, 0, 0); // just to force alpha array allocation with
                                  // the right size
        Arrays.fill(imData.alphaData, (byte) 0); // set whole image as
                                                 // transparent

        Image image = new Image(display, imData);

        Color fill = config.getDCElementExternalFillColor();
        GC gc = new GC(image);
        gc.setBackground(fill);
        gc.setAlpha(alpha);
        fillRectangle(gc, scale, 4, 0, 2, 2);
        fillRectangle(gc, scale, 2, 2, 2, 2);
        fillRectangle(gc, scale, 0, 4, 2, 2);
        // CHECKSTYLEON MagicNumber

        gc.dispose();

        return new Pattern(display, image);
    }

    /**
     * Get a pattern for the extensions in polyominoes.
     * 
     * @param alpha
     *            transparency
     * @return the pattern
     */
    public Pattern getPolyominoExtensionPattern(final int alpha) {
        // This is hacky:
        // http://stackoverflow.com/questions/8270474/drawing-on-a-transparent-image-using-java-swt
        // but we want some totally transparent pixels
        // allocate an image data
        // CHECKSTYLEOFF MagicNumber
        ImageData imData = new ImageData(6 * scale, 6 * scale, 24, new PaletteData(0xff0000, 0x00ff00, 0x0000ff));
        imData.setAlpha(0, 0, 0); // just to force alpha array allocation with
                                  // the right size
        Arrays.fill(imData.alphaData, (byte) 0); // set whole image as
                                                 // transparent

        Image image = new Image(display, imData);

        Color fill = config.getPolyominoWeaklyBlockedColor();
        GC gc = new GC(image);
        gc.setBackground(fill);
        gc.setAlpha(alpha);
        fillRectangle(gc, scale, 0, 0, 2, 2);
        fillRectangle(gc, scale, 2, 2, 2, 2);
        fillRectangle(gc, scale, 4, 4, 2, 2);
        // CHECKSTYLEON MagicNumber

        gc.dispose();

        if (scale != 1) {
            imData = image.getImageData();
            image.dispose();
            imData.scaledTo(imData.width * scale, imData.height * scale);
            image = new Image(display, imData);
        }

        return new Pattern(display, image);
    }

    /**
     * Get a pattern for highlighting the center of a polyomino.
     * 
     * @param alpha
     *            transparency
     * @return the pattern
     */
    public Pattern getPolyominoCenterPattern(final int alpha) {
        // This is hacky:
        // http://stackoverflow.com/questions/8270474/drawing-on-a-transparent-image-using-java-swt
        // but we want some totally transparent pixels
        // allocate an image data
        // CHECKSTYLEOFF MagicNumber
        ImageData imData = new ImageData(6 * scale, 6 * scale, 24, new PaletteData(0xff0000, 0x00ff00, 0x0000ff));
        imData.setAlpha(0, 0, 0); // just to force alpha array allocation with
                                  // the right size
        Arrays.fill(imData.alphaData, (byte) 0); // set whole image as
                                                 // transparent

        Image image = new Image(display, imData);

        Color fill = config.getPolyominoCenterColor();
        GC gc = new GC(image);
        gc.setBackground(fill);
        gc.setAlpha(alpha);
        fillRectangle(gc, scale, 2, 0, 2, 2);
        fillRectangle(gc, scale, 0, 2, 2, 2);
        fillRectangle(gc, scale, 4, 2, 2, 2);
        fillRectangle(gc, scale, 2, 4, 2, 2);
        // CHECKSTYLEON MagicNumber

        gc.dispose();

        return new Pattern(display, image);
    }

    /**
     * Fills a rectangle inside the pattern canvas.
     * 
     * @param gc
     *            Graphical context
     * @param scaling
     *            Scaling factor
     * @param x
     *            X-Coordinate
     * @param y
     *            Y-Coordinate
     * @param width
     *            With of area to be colored
     * @param height
     *            Height of area to be colored
     */
    private void fillRectangle(final GC gc, final int scaling, final int x, final int y, final int width,
            final int height) {
        gc.fillRectangle(x * scaling, y * scaling, width * scaling, height * scaling);
    }
}
