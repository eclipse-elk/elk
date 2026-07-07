/*******************************************************************************
 * Copyright (c) 2018, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.rectpacking.p1widthapproximation;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.rectpacking.options.OptimizationGoal;
import org.eclipse.elk.alg.rectpacking.util.DrawingData;
import org.eclipse.elk.alg.rectpacking.util.DrawingDataDescriptor;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.graph.ElkNode;

import com.google.common.collect.Lists;

/**
 * Class that handles the first iteration of the algorithm, producing an approximated bounding box for the packing.
 */
public class AreaApproximation {
    //////////////////////////////////////////////////////////////////
    // Fields.
    /** Desired aspect ratio. */
    private double aspectRatio;
    /** Optimization goal. */
    private OptimizationGoal goal;
    /** Shift when placing behind or below the last placed rectangle. */
    private boolean lpShift;

    //////////////////////////////////////////////////////////////////
    // Constructors.
    /**
     * Constructs an object that the first iteration can be executed on.
     * 
     * @param aspectRatio The desired aspect ratio.
     * @param goal The optimization goal.
     * @param lpShift Whether a shift should happen after placing near the last placed node.
     */
    public AreaApproximation(final double aspectRatio, final OptimizationGoal goal, final boolean lpShift) {
        this.aspectRatio = aspectRatio;
        this.goal = goal;
        this.lpShift = lpShift;
    }


    //////////////////////////////////////////////////////////////////
    // Public methods.
    /**
     * Calculates a drawing for the given rectangles according options set in the object creation of
     * {@link AreaApproximation}. This method also sets the coordinates for the rectangles.
     * 
     * @param rectangles The rectangles to place.
     * @param nodeNodeSpacing The spacing between two nodes.
     * @return A drawing calculated by this methods algorithm.
     */
    public DrawingData approxBoundingBox(final List<ElkNode> rectangles, final double nodeNodeSpacing,
            final ElkPadding padding) {
        // Place first box.
        ElkNode firstRect = rectangles.get(0);
        firstRect.setX(0);
        firstRect.setY(0);
        List<ElkNode> placedRects = new ArrayList<>();
        placedRects.add(firstRect);
        ElkNode lastPlaced = firstRect;
        DrawingData currentValues = new DrawingData(this.aspectRatio, firstRect.getWidth(), firstRect.getHeight(),
                DrawingDataDescriptor.WHOLE_DRAWING);

        // Place the other boxes.
        for (int rectangleIdx = 1; rectangleIdx < rectangles.size(); rectangleIdx++) {
            ElkNode toPlace = rectangles.get(rectangleIdx);

            // Determine drawing metrics for different candidate positions/placement options
            DrawingData opt1 = calcValuesForOpt(DrawingDataDescriptor.CANDIDATE_POSITION_LAST_PLACED_RIGHT, toPlace,
                    lastPlaced, currentValues, placedRects, nodeNodeSpacing);

            DrawingData opt2 = calcValuesForOpt(DrawingDataDescriptor.CANDIDATE_POSITION_LAST_PLACED_BELOW, toPlace,
                    lastPlaced, currentValues, placedRects, nodeNodeSpacing);

            DrawingData opt3 = calcValuesForOpt(DrawingDataDescriptor.CANDIDATE_POSITION_WHOLE_DRAWING_RIGHT, toPlace,
                    lastPlaced, currentValues, placedRects, nodeNodeSpacing);

            DrawingData opt4 = calcValuesForOpt(DrawingDataDescriptor.CANDIDATE_POSITION_WHOLE_DRAWING_BELOW, toPlace,
                    lastPlaced, currentValues, placedRects, nodeNodeSpacing);

            DrawingData bestOpt = findBestCandidate(opt1, opt2, opt3, opt4, toPlace, lastPlaced, padding);

            toPlace.setX(bestOpt.getNextXcoordinate());
            toPlace.setY(bestOpt.getNextYcoordinate());
            bestOpt.setPlacementOption(DrawingDataDescriptor.WHOLE_DRAWING);
            currentValues = bestOpt;
            lastPlaced = toPlace;
            placedRects.add(toPlace);
        }

        return currentValues;
    }

    //////////////////////////////////////////////////////////////////
    // Helper methods.

    /**
     * Determines according to the selected optimization goal, which option out of the given ones is the best.
     * 
     * @param opt1 The option lastPlaced right.
     * @param opt2 The option lastPlaced bottom.
     * @param opt3 The option complete drawing right.
     * @param opt4 The option complete drawing below.
     * @param toPlace The rectangle to be placed on drawing area.
     * @param lastPlaced The rectangle that was placed last on the drawing area.
     * @return Returns the best option out of the given ones or null, if there is no best option found.
     */
    private DrawingData findBestCandidate(final DrawingData opt1, final DrawingData opt2, final DrawingData opt3,
            final DrawingData opt4, final ElkNode toPlace, final ElkNode lastPlaced, final ElkPadding padding) {
        List<DrawingData> candidates = Lists.newArrayList(opt1, opt2, opt3, opt4);
        List<BestCandidateFilter> filters = null;

        switch (goal) {
        // Sets the order of the filters according to the given goal.
        case MAX_SCALE_DRIVEN:
            filters = Lists.newArrayList(new ScaleMeasureFilter(), new AreaFilter(), new AspectRatioFilter());
            break;
        case ASPECT_RATIO_DRIVEN:
            filters = Lists.newArrayList(new AspectRatioFilter(), new AreaFilter(), new ScaleMeasureFilter());
            break;
        case AREA_DRIVEN:
            filters = Lists.newArrayList(new AreaFilter(), new ScaleMeasureFilter(), new AspectRatioFilter());
            break;
        default:
        }

        // Filter the candidates according to the order of the filters using the strategy pattern.
        for (BestCandidateFilter filter : filters) {
            if (candidates.size() > 1) {
                candidates = filter.filterList(candidates, aspectRatio, padding);
            }
        }

        // Only one candidate remains.
        if (candidates.size() == 1) {
            return candidates.get(candidates.size() - 1);
        }
        // Multiple options have the same value for every benchmark. These special cases are caught in the following.
        if (candidates.size() == 2) {
            return checkSpecialCases(candidates.get(0), candidates.get(1), lastPlaced, toPlace);
        }
        return null;
    }

    /**
     * Checks the two options and their possible combinations after being sorted out be the packing strategies.
     * Determines the best option depending on the options given.
     * 
     * @param drawing1 A drawing to compare.
     * @param drawing2 A drawing to compare.
     * @param lastPlaced The rectangle that was last placed before the two drawings were calculated.
     * @param toPlace The rectangle that was newly placed in the two drawings.
     * @return The better drawing out of the two.
     */
    private DrawingData checkSpecialCases(final DrawingData drawing1, final DrawingData drawing2,
            final ElkNode lastPlaced, final ElkNode toPlace) {
        DrawingDataDescriptor firstOpt = drawing1.getPlacementOption();
        DrawingDataDescriptor secondOpt = drawing2.getPlacementOption();

        boolean firstOptLPBorCDB = (firstOpt == DrawingDataDescriptor.CANDIDATE_POSITION_LAST_PLACED_BELOW
                || firstOpt == DrawingDataDescriptor.CANDIDATE_POSITION_WHOLE_DRAWING_BELOW);
        boolean secondOptLPBorCDB = (secondOpt == DrawingDataDescriptor.CANDIDATE_POSITION_LAST_PLACED_BELOW
                || secondOpt == DrawingDataDescriptor.CANDIDATE_POSITION_WHOLE_DRAWING_BELOW);

        boolean firstOptLPRorCDR = (firstOpt == DrawingDataDescriptor.CANDIDATE_POSITION_LAST_PLACED_RIGHT
                || firstOpt == DrawingDataDescriptor.CANDIDATE_POSITION_WHOLE_DRAWING_RIGHT);
        boolean secondOptLPRorCDR = (secondOpt == DrawingDataDescriptor.CANDIDATE_POSITION_LAST_PLACED_RIGHT
                || secondOpt == DrawingDataDescriptor.CANDIDATE_POSITION_WHOLE_DRAWING_RIGHT);

        boolean firstOptLPRorLPB = (firstOpt == DrawingDataDescriptor.CANDIDATE_POSITION_LAST_PLACED_RIGHT
                || firstOpt == DrawingDataDescriptor.CANDIDATE_POSITION_LAST_PLACED_BELOW);
        boolean secondOptLPRorLPB = (secondOpt == DrawingDataDescriptor.CANDIDATE_POSITION_LAST_PLACED_RIGHT
                || secondOpt == DrawingDataDescriptor.CANDIDATE_POSITION_LAST_PLACED_BELOW);

        if (firstOptLPBorCDB && secondOptLPBorCDB) {
            // If placing it LPB and WDB produces the same values. Take WDB.
            if (drawing1.getPlacementOption() == DrawingDataDescriptor.CANDIDATE_POSITION_WHOLE_DRAWING_BELOW) {
                return drawing1;
            } else {
                return drawing2;
            }
        } else if (firstOptLPRorCDR && secondOptLPRorCDR) {
            // If placing it LPR and WDR produces the same values. Take WDR.
            if (drawing1.getPlacementOption() == DrawingDataDescriptor.CANDIDATE_POSITION_WHOLE_DRAWING_RIGHT) {
                return drawing1;
            } else {
                return drawing2;
            }
        } else if (firstOptLPRorLPB && secondOptLPRorLPB) {
            // If LPR AND LPB produce the same values. Take the option producing less area with the last placed
            // rectangle and rectangle to place.
            DrawingData lprOpt;
            DrawingData lpbOpt;
            if (firstOpt == DrawingDataDescriptor.CANDIDATE_POSITION_LAST_PLACED_RIGHT) {
                lprOpt = drawing1;
                lpbOpt = drawing2;
            } else {
                lprOpt = drawing2;
                lpbOpt = drawing1;
            }

            double areaLPR = Calculations.calculateAreaLPR(lastPlaced, toPlace, lprOpt);
            double areaLPB = Calculations.calculateAreaLPB(lastPlaced, toPlace, lpbOpt);

            if (areaLPR <= areaLPB) {
                if (drawing1.getPlacementOption() == DrawingDataDescriptor.CANDIDATE_POSITION_LAST_PLACED_RIGHT) {
                    return drawing1;
                } else {
                    return drawing2;
                }
            } else {
                if (drawing1.getPlacementOption() == DrawingDataDescriptor.CANDIDATE_POSITION_LAST_PLACED_BELOW) {
                    return drawing1;
                } else {
                    return drawing2;
                }
            }
        }

        return drawing1;
    }

    /**
     * Calculates drawing data for the given parameters including x and y coordinate for the rectangle toPlace.
     * 
     * @param option The placement option the values are calculated for.
     * @param toPlace The rectangle to be placed.
     * @param lastPlaced The rectangle that was placed last.
     * @param drawing The current drawing containing width and height, besides others.
     * @param placedRects The already placed rectangles.
     * @param nodeNodeSpacing The spacing between two nodes.
     * @return An {@link DrawingData} object containing the values after the rectangle would be placed.
     */
    private DrawingData calcValuesForOpt(final DrawingDataDescriptor option, final ElkNode toPlace,
            final ElkNode lastPlaced, final DrawingData drawing, final List<ElkNode> placedRects,
            final double nodeNodeSpacing) {

        double x = 0;
        double y = 0;
        double drawingWidth = drawing.getDrawingWidth();
        double drawingHeight = drawing.getDrawingHeight();
        double heightToPlace = toPlace.getHeight();
        double widthToPlace = toPlace.getWidth();

        double width, height;
        switch (option) {
        case CANDIDATE_POSITION_LAST_PLACED_RIGHT:
            x = lastPlaced.getX() + lastPlaced.getWidth() + nodeNodeSpacing;
            if (lpShift) {
                y = Calculations.calculateYforLPR(x, placedRects, lastPlaced, nodeNodeSpacing);
            } else {
                y = lastPlaced.getY();
            }

            width = Calculations.getWidthLPRorLPB(drawingWidth, x, widthToPlace);
            height = Calculations.getHeightLPRorLPB(drawingHeight, y, heightToPlace);
            break;
            
        case CANDIDATE_POSITION_LAST_PLACED_BELOW:
            y = lastPlaced.getY() + lastPlaced.getHeight() + nodeNodeSpacing;
            if (lpShift) {
                x = Calculations.calculateXforLPB(y, placedRects, lastPlaced, nodeNodeSpacing);
            } else {
                x = lastPlaced.getX();
            }

            width = Calculations.getWidthLPRorLPB(drawingWidth, x, widthToPlace);
            height = Calculations.getHeightLPRorLPB(drawingHeight, y, heightToPlace);
            break;
            
        case CANDIDATE_POSITION_WHOLE_DRAWING_RIGHT:
            x = drawingWidth + nodeNodeSpacing;
            y = 0;

            width = drawingWidth + nodeNodeSpacing + widthToPlace;
            height = Math.max(drawingHeight, heightToPlace);
            break;
            
        case CANDIDATE_POSITION_WHOLE_DRAWING_BELOW:
            x = 0;
            y = drawingHeight + nodeNodeSpacing;

            width = Math.max(drawingWidth, widthToPlace);
            height = drawingHeight + nodeNodeSpacing + heightToPlace;
            break;
        default:
            throw new IllegalArgumentException("IllegalPlacementOption.");
        }

        // LPR and LPB have more values to be set in newDrawing.
        DrawingData newDrawing = new DrawingData(aspectRatio, width, height, option, x, y);
        return newDrawing;
    }
}
