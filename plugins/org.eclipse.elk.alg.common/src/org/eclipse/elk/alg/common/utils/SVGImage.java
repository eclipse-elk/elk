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
package org.eclipse.elk.alg.common.utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;
import java.util.Map.Entry;

import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;

import com.google.common.collect.Maps;

/**
 * SVG output for debugging purposes.
 * 
 * Note: when adding/removing methods from this class, make sure to adjust the emulating class in elkjs.
 */
public class SVGImage {
    // CHECKSTYLEOFF MagicNumber
    private final double spacing = 20;
    private double viewBoxX = Double.POSITIVE_INFINITY, viewBoxY = Double.POSITIVE_INFINITY, viewBoxW = 0, viewBoxH = 0,
            viewBoxX2 = Double.NEGATIVE_INFINITY, viewBoxY2 = Double.NEGATIVE_INFINITY;
    private boolean enableUpdateViewBox = true;
    private String current = "main";
    private Map<String, String> groups = Maps.newLinkedHashMap();
    // SUPPRESS CHECKSTYLE NEXT 8 VisibilityModifier
    /** Path to the file where the image should be saved. The file extension
     *  and a possible numeration are appended automatically. */
    public String fileName = null;
    private int index = 0;
    /** Enables the outpout. */
    public boolean debug = false;
    /** The image will be scaled down if it exceeds this size. */
    public double maxArea = 8294400; // 4k

    /**
     * The constructor disables the output if the filename is null.
     * 
     * @param file the SVG's filename
     */
    public SVGImage(final String file) {
        if (file != null && file.length() > 0 && !file.startsWith("null")) {
            fileName = file;
            debug = true;
        }
        groups.put(current, "");
    }
    
    /**
     * Select an SVG group.
     * 
     * @param key the identifier of an existing or a new group
     * @return the image with the selected group being active
     */
    public SVGImage g(final String key) {
        if (debug) {
            if (!groups.containsKey(key)) {
                groups.put(key, "");
            }
            current = key;
        }
        return this;
    }

    /**
     * Add SVG groups.
     * 
     * @param keys the new groups' identifiers
     */
    public void addGroups(final String... keys) {
        if (debug) {
            for (String key : keys) {
                groups.put(key, ""); // clears existing
            }
        }
    }
    
    /**
     * Clears the specified group.
     * 
     * @param key identifier of the group
     */
    public void clearGroup(final String key) {
        if (debug) {
            groups.replace(key, "");
        }
    }

    /**
     * Removes group.
     * 
     * @param key identifier of the group
     */
    public void removeGroup(final String key) {
        if (debug) {
            if (key.equals("main")) {
                clearGroup(key);
            } else {
                groups.remove(key);
            }
        }
    }

    /**
     * Sets the SVG's viewbox to a fixed size and position and disables its
     * automatic updating.
     * 
     * @param x x
     * @param y y
     * @param w width
     * @param h height
     */
    public void setViewBox(final double x, final double y, final double w, final double h) {
        if (debug) {
            viewBoxX = x;
            viewBoxY = y;
            viewBoxW = w;
            viewBoxH = h;
            enableUpdateViewBox = false;
        }
    }

    /**
     * Clears the entire image and resets the output according to the 
     * availability of a filename.
     */
    public void clear() {
        debug = fileName != null && fileName.length() > 0;
        if (debug) {
            groups.clear();
            current = "main";
            groups.put(current, "");
            if (enableUpdateViewBox) {
                viewBoxX = Double.POSITIVE_INFINITY;
                viewBoxY = Double.POSITIVE_INFINITY;
                viewBoxW = 0;
                viewBoxH = 0;
                viewBoxX2 = Double.NEGATIVE_INFINITY;
                viewBoxY2 = Double.NEGATIVE_INFINITY;
            }
        }
    }
    
    // CHECKSTYLEOFF Javadoc

    public void addElementStr(final String element) {
        if (debug) {
            groups.replace(current, groups.get(current).concat(element + "\n"));
            current = "main";
        }
    }

    public void addCircle(final double x, final double y) {
        if (debug) {
            addCircle(x, y, 5, "stroke=\"black\" stroke-width=\"1\" fill=\"none\"");
        }
    }

    public void addCircle(final double x, final double y, final double r, final String attributes) {
        if (debug) {
            addElementStr("<circle cx=\"" + x + "\" cy=\"" + y + "\" r=\"" + r + "\" " + attributes + " />");
            updateViewBox(new KVector(x - r, y - r), new KVector(x + r, y + r));
        }
    }

    public void addLine(final double x1, final double y1, final double x2, final double y2) {
        if (debug) {
            addLine(x1, y1, x2, y2, "stroke=\"black\" stroke-width=\"1\"");
        }
    }

    public void addLine(final double x1, final double y1, final double x2, final double y2, final String attributes) {
        if (debug) {
            addElementStr("<line x1=\"" + x1 + "\" y1=\"" + y1 + "\" x2=\"" + x2 + "\" y2=\"" + y2 + "\" " + attributes
                    + " />");
            updateViewBox(new KVector(x1, y1), new KVector(x2, y2));
        }
    }

    public void addRect(final double x, final double y, final double w, final double h, final String attributes) {
        if (debug) {
            addElementStr("<rect x=\"" + x + "\" y=\"" + y + "\" width=\"" + w + "\" height=\"" + h + "\" " + attributes
                    + " />");
            updateViewBox(new KVector(x, y), new KVector(x + w, y + h));
        }
    }

    public void addRect(final ElkRectangle r, final String attributes) {
        if (debug) {
            addRect(r.x, r.y, r.width, r.height, attributes);
        }
    }

    public void addPoly(final String attributes, final KVector... points) {
        if (debug) {
            String str = "<polyline points=\"";
            for (KVector p : points) {
                str += p.x + "," + p.y + " ";
            }
            addElementStr(str + "\" " + attributes + " />");
            updateViewBox(points);
        }
    }
    
    public void addText(final double x, final double y, final String text, final String attributes) {
    // anchor bottom left
        if (debug) {
            addElementStr("<text x=\"" + x + "\" y=\"" + y + "\" " + attributes + ">"
                    + text + "</text>");
            // no update because size unknown
        }
    }

    private void updateViewBox(final KVector... points) {
        if (debug) {
            if (enableUpdateViewBox) {
                for (KVector point : points) {
                    if (Double.isFinite(point.x) && Double.isFinite(point.y)) {
                        viewBoxX = Math.min(viewBoxX, point.x);
                        viewBoxY = Math.min(viewBoxY, point.y);
                        viewBoxX2 = Math.max(viewBoxX2, point.x);
                        viewBoxY2 = Math.max(viewBoxY2, point.y);
                    }
                }
                viewBoxW = viewBoxX2 - viewBoxX;
                viewBoxH = viewBoxY2 - viewBoxY;
            }
        }
    }

    /**
     * Writes the SVG to the given file.
     * 
     * @param fName filename without file extension
     */
    public void save(final String fName) {
        if (debug) {
            PrintWriter out;
            try {
                double area = viewBoxH * viewBoxW;
                double scale = 1;
                if (area > maxArea) {
                    scale = Math.sqrt(maxArea / area);
                }
                new File(fName).getParentFile().mkdirs();
                out = new PrintWriter(new FileWriter(fName + ".svg"));

                out.println("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
                out.println("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\"" + "  viewBox=\""
                        + (viewBoxX * scale - spacing) + " " + (viewBoxY * scale - spacing) + " "
                        + (viewBoxW * scale + 2 * spacing) + " " + (viewBoxH * scale + 2 * spacing) + "\">");
                for (Entry<String, String> entry : groups.entrySet()) {
                    out.println("<g transform=\"scale(" + scale + ")\" id=\"" + entry.getKey() + "\">\n" 
                                + entry.getValue() + "</g>");
                }
                out.println("</svg>");

                out.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Writes the SVG.
     */
    public void save() {
        save(fileName);
    }

    /**
     * Writes the SVG to the given file and appends an incrementing number.
     * 
     * @param fName filename without file extension
     */
    public void isave(final String fName) {
        save(fName + String.format("%1$03d", index));
        index++;
    }

    /**
     * Writes the SVG and appends an incrementing number.
     */
    public void isave() {
        isave(fileName);
    }
}
