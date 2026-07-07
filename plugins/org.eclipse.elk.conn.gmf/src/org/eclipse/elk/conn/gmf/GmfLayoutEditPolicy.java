/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.conn.gmf;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.eclipse.draw2d.ConnectionAnchor;
import org.eclipse.draw2d.ConnectionLocator;
import org.eclipse.draw2d.IFigure;
import org.eclipse.draw2d.geometry.Dimension;
import org.eclipse.draw2d.geometry.Point;
import org.eclipse.draw2d.geometry.PointList;
import org.eclipse.draw2d.geometry.PrecisionPoint;
import org.eclipse.draw2d.geometry.Rectangle;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.EdgeLabelPlacement;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.core.util.WrappedException;
import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.ElkShape;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.ecore.EReference;
import org.eclipse.gef.GraphicalEditPart;
import org.eclipse.gef.Request;
import org.eclipse.gef.commands.Command;
import org.eclipse.gef.editpolicies.AbstractEditPolicy;
import org.eclipse.gmf.runtime.diagram.core.util.ViewUtil;
import org.eclipse.gmf.runtime.diagram.ui.commands.ICommandProxy;
import org.eclipse.gmf.runtime.diagram.ui.editparts.ConnectionEditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.IGraphicalEditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.INodeEditPart;
import org.eclipse.gmf.runtime.diagram.ui.editparts.LabelEditPart;
import org.eclipse.gmf.runtime.draw2d.ui.geometry.LineSeg;
import org.eclipse.gmf.runtime.draw2d.ui.geometry.PointListUtilities;
import org.eclipse.gmf.runtime.emf.core.util.EObjectAdapter;
import org.eclipse.gmf.runtime.gef.ui.figures.SlidableAnchor;
import org.eclipse.gmf.runtime.notation.Edge;
import org.eclipse.gmf.runtime.notation.NotationPackage;
import org.eclipse.gmf.runtime.notation.View;

/**
 * Edit policy used to apply layout. This edit policy creates a {@link GmfLayoutCommand} to directly
 * manipulate layout data in the GMF notation model.
 * 
 * @author msp
 * @see org.eclipse.gmf.runtime.diagram.ui.editpolicies.XYLayoutEditPolicy
 * @see org.eclipse.gmf.runtime.diagram.ui.editpolicies.ConnectionBendpointEditPolicy
 * @see org.eclipse.gmf.runtime.diagram.ui.editpolicies.GraphicalNodeEditPolicy
 */
public class GmfLayoutEditPolicy extends AbstractEditPolicy {

    /** map of edge layouts to existing point lists. */
    private Map<ElkEdgeSection, PointList> pointListMap = new HashMap<>();

    @Override
    public boolean understandsRequest(final Request req) {
        return (ApplyLayoutRequest.REQ_APPLY_LAYOUT.equals(req.getType()));
    }

    @Override
    public Command getCommand(final Request request) {
        if (ApplyLayoutRequest.REQ_APPLY_LAYOUT.equals(request.getType())) {
            if (request instanceof ApplyLayoutRequest) {
                ApplyLayoutRequest layoutRequest = (ApplyLayoutRequest) request;
                IGraphicalEditPart hostEditPart = (IGraphicalEditPart) getHost();
                GmfLayoutCommand command = new GmfLayoutCommand(hostEditPart.getEditingDomain(),
                        "Automatic Layout", new EObjectAdapter((View) hostEditPart.getModel()));
                double scale = layoutRequest.getScale();

                // retrieve layout data from the request and compute layout data for the command
                for (Pair<ElkGraphElement, GraphicalEditPart> layoutPair : layoutRequest.getElements()) {
                    if (layoutPair.getFirst() instanceof ElkNode) {
                        addShapeLayout(command, (ElkShape) layoutPair.getFirst(), layoutPair.getSecond(), scale);
                    } else if (layoutPair.getFirst() instanceof ElkPort) {
                        addShapeLayout(command, (ElkPort) layoutPair.getFirst(), layoutPair.getSecond(), scale);
                    } else if (layoutPair.getFirst() instanceof ElkEdge) {
                        addEdgeLayout(command, (ElkEdge) layoutPair.getFirst(),
                                (ConnectionEditPart) layoutPair.getSecond(), scale);
                    } else if (layoutPair.getFirst() instanceof ElkLabel) {
                        addLabelLayout(command, (ElkLabel) layoutPair.getFirst(), layoutPair.getSecond(), scale);
                    }
                }

                // TODO Make this configurable?
                command.setObliqueRouting(true);

                pointListMap.clear();
                return new ICommandProxy(command);
            } else {
                return null;
            }
        } else {
            return super.getCommand(request);
        }
    }

    /**
     * Adds a shape layout to the given command.
     * 
     * @param command
     *            command to which a shape layout shall be added
     * @param elkShape
     *            graph element with layout data
     * @param editPart
     *            edit part to which layout is applied
     * @param scale
     *            scale factor for coordinates
     */
    private void addShapeLayout(final GmfLayoutCommand command, final ElkShape elkShape,
            final GraphicalEditPart editPart, final double scale) {
        
        View view = (View) editPart.getModel();

        // check whether the location has changed
        Point newLocation = new Point((int) (elkShape.getX() * scale), (int) (elkShape.getY() * scale));
        Object oldx = ViewUtil.getStructuralFeatureValue(view, NotationPackage.eINSTANCE.getLocation_X());
        Object oldy = ViewUtil.getStructuralFeatureValue(view, NotationPackage.eINSTANCE.getLocation_Y());
        
        if (oldx != null && oldy != null && newLocation.x == (Integer) oldx && newLocation.y == (Integer) oldy) {
            newLocation = null;
        }

        // check whether the size has changed
        Dimension newSize = new Dimension(
                (int) (elkShape.getWidth() * scale),
                (int) (elkShape.getHeight() * scale));
        Object oldWidth = ViewUtil.getStructuralFeatureValue(view, NotationPackage.eINSTANCE.getSize_Width());
        Object oldHeight = ViewUtil.getStructuralFeatureValue(view, NotationPackage.eINSTANCE.getSize_Height());
        
        if (oldWidth != null && oldHeight != null && newSize.width == (Integer) oldWidth
                && newSize.height == (Integer) oldHeight) {
            
            newSize = null;
        }

        if (newLocation != null || newSize != null) {
            command.addShapeLayout(view, newLocation, newSize);
        }
    }

    /**
     * Adds an edge layout to the given command.
     * 
     * @param command
     *            command to which an edge layout shall be added
     * @param elkEdge
     *            edge with layout data
     * @param connectionEditPart
     *            edit part to which layout is applied
     * @param scale
     *            scale factor for coordinates
     */
    private void addEdgeLayout(final GmfLayoutCommand command, final ElkEdge elkEdge,
            final ConnectionEditPart connectionEditPart, final double scale) {
        
        if (connectionEditPart.getSource() != null && connectionEditPart.getTarget() != null) {
            // create source terminal identifier
            INodeEditPart sourceEditPart = (INodeEditPart) connectionEditPart.getSource();
            ConnectionAnchor sourceAnchor;
            if (sourceEditPart instanceof ConnectionEditPart) {
                // if the edge source is a connection, don't consider the source point
                sourceAnchor = new SlidableAnchor(sourceEditPart.getFigure());
            } else {
                KVector sourceRel = getRelativeSourcePoint(elkEdge);
                sourceAnchor = new SlidableAnchor(sourceEditPart.getFigure(),
                        new PrecisionPoint(sourceRel.x, sourceRel.y));
            }
            String sourceTerminal = sourceEditPart.mapConnectionAnchorToTerminal(sourceAnchor);

            // create target terminal identifier
            INodeEditPart targetEditPart = (INodeEditPart) connectionEditPart.getTarget();
            ConnectionAnchor targetAnchor;
            if (targetEditPart instanceof ConnectionEditPart) {
                // if the edge target is a connection, don't consider the target point
                targetAnchor = new SlidableAnchor(targetEditPart.getFigure());
            } else {
                KVector targetRel = getRelativeTargetPoint(elkEdge);
                targetAnchor = new SlidableAnchor(targetEditPart.getFigure(),
                        new PrecisionPoint(targetRel.x, targetRel.y));
            }
            String targetTerminal = targetEditPart.mapConnectionAnchorToTerminal(targetAnchor);

            PointList bendPoints = getBendPoints(elkEdge, connectionEditPart.getFigure(), scale);

            // check whether the connection is a note attachment to an edge, then remove bend points
            if (sourceEditPart instanceof ConnectionEditPart
                    || targetEditPart instanceof ConnectionEditPart) {
                while (bendPoints.size() > 2) {
                    bendPoints.removePoint(1);
                }
            }
            
            // retrieve junction points and transform them to absolute coordinates
            KVectorChain junctionPoints = elkEdge.getProperty(CoreOptions.JUNCTION_POINTS);
            String serializedJP = null;
            if (junctionPoints != null) {
                for (KVector point : junctionPoints) {
                    ElkUtil.toAbsolute(point, elkEdge.getContainingNode());
                }
                serializedJP = junctionPoints.toString();
            }

            command.addEdgeLayout((Edge) connectionEditPart.getModel(), bendPoints, sourceTerminal,
                    targetTerminal, serializedJP);
        }
    }

    /**
     * Create a vector that contains the relative position of the source point to the corresponding source node or
     * port.
     * 
     * @param edge
     *            an edge
     * @return the relative source point
     */
    private KVector getRelativeSourcePoint(final ElkEdge edge) {
        // The edge should have exactly one source shape
        ElkConnectableShape sourceShape = edge.getSources().get(0);
        
        // The edge should have one edge section after layout
        ElkEdgeSection edgeSection = edge.getSections().get(0);
        KVector sourcePoint = new KVector(edgeSection.getStartX(), edgeSection.getStartY());
        
        // We will now make the source point absolute, and then relative to the source node
        ElkUtil.toAbsolute(sourcePoint, edge.getContainingNode());
        ElkUtil.toRelative(sourcePoint, ElkGraphUtil.connectableShapeToNode(sourceShape));
        
        // The end result will be coordinates between 0 and 1, with 0 being at the left / top or the source shape and
        // 1 being at the right / bottom
        if (sourceShape instanceof ElkPort) {
            ElkPort sourcePort = (ElkPort) sourceShape;
            
            // calculate the relative position to the port size
            if (sourcePort.getWidth() <= 0) {
                sourcePoint.x = 0;
            } else {
                sourcePoint.x = (sourcePoint.x - sourcePort.getX()) / sourcePort.getWidth();
            }
            
            if (sourcePort.getHeight() <= 0) {
                sourcePoint.y = 0;
            } else {
                sourcePoint.y = (sourcePoint.y - sourcePort.getY()) / sourcePort.getHeight();
            }
        } else {
            // calculate the relative position to the node size
            if (sourceShape.getWidth() <= 0) {
                sourcePoint.x = 0;
            } else {
                sourcePoint.x /= sourceShape.getWidth();
            }
            if (sourceShape.getHeight() <= 0) {
                sourcePoint.y = 0;
            } else {
                sourcePoint.y /= sourceShape.getHeight();
            }
        }

        // check the bound of the relative position
        return sourcePoint.bound(0, 0, 1, 1);
    }

    /**
     * Create a vector that contains the relative position of the target point to the corresponding
     * target node or port.
     * 
     * @param edge
     *            an edge
     * @return the relative target point
     */
    private KVector getRelativeTargetPoint(final ElkEdge edge) {
        // The edge should have exactly one source shape
        ElkConnectableShape targetShape = edge.getTargets().get(0);
        
        // The edge should have one edge section after layout
        ElkEdgeSection edgeSection = edge.getSections().get(0);
        KVector targetPoint = new KVector(edgeSection.getEndX(), edgeSection.getEndY());
        
        // We will now make the source point absolute, and then relative to the source node
        ElkUtil.toAbsolute(targetPoint, edge.getContainingNode());
        ElkUtil.toRelative(targetPoint, ElkGraphUtil.connectableShapeToNode(targetShape));
        
        // The end result will be coordinates between 0 and 1, with 0 being at the left / top or the source shape and
        // 1 being at the right / bottom
        if (targetShape instanceof ElkPort) {
            ElkPort targetPort = (ElkPort) targetShape;
            
            // calculate the relative position to the port size
            if (targetPort.getWidth() <= 0) {
                targetPoint.x = 0;
            } else {
                targetPoint.x = (targetPoint.x - targetPort.getX()) / targetPort.getWidth();
            }
            
            if (targetPort.getHeight() <= 0) {
                targetPoint.y = 0;
            } else {
                targetPoint.y = (targetPoint.y - targetPort.getY()) / targetPort.getHeight();
            }
        } else {
            // calculate the relative position to the node size
            if (targetShape.getWidth() <= 0) {
                targetPoint.x = 0;
            } else {
                targetPoint.x /= targetShape.getWidth();
            }
            if (targetShape.getHeight() <= 0) {
                targetPoint.y = 0;
            } else {
                targetPoint.y /= targetShape.getHeight();
            }
        }

        // check the bound of the relative position
        return targetPoint.bound(0, 0, 1, 1);
    }

    /** see LabelViewConstants.TARGET_LOCATION. */
    private static final int SOURCE_LOCATION = 85;
    /** see LabelViewConstants.MIDDLE_LOCATION. */
    private static final int MIDDLE_LOCATION = 50;
    /** see LabelViewConstants.SOURCE_LOCATION. */
    private static final int TARGET_LOCATION = 15;

    /**
     * Adds an edge label layout to the given command.
     * 
     * @param command
     *            command to which the edge label layout shall be added
     * @param klabel
     *            label with layout data
     * @param labelEditPart
     *            edit part to which layout is applied
     * @param scale
     *            scale factor for coordinates
     */
    private void addLabelLayout(final GmfLayoutCommand command, final ElkLabel klabel,
            final GraphicalEditPart labelEditPart, final double scale) {
        
        ElkGraphElement parent = klabel.getParent();
        
        // node and port labels are processed separately
        if (parent instanceof ElkNode || parent instanceof ElkPort) {
            View view = (View) labelEditPart.getModel();
            int xpos = (int) (klabel.getX() * scale);
            int ypos = (int) (klabel.getY() * scale);
            Object oldx = ViewUtil.getStructuralFeatureValue(view, NotationPackage.eINSTANCE.getLocation_X());
            Object oldy = ViewUtil.getStructuralFeatureValue(view, NotationPackage.eINSTANCE.getLocation_Y());
            
            if (oldx == null || oldy == null || xpos != (Integer) oldx || ypos != (Integer) oldy) {
                command.addShapeLayout(view, new Point(xpos, ypos), null);
            }
            return;
        } else if (parent instanceof ElkEdge) {
            // calculate direct new location of the label
            Rectangle targetBounds = new Rectangle(labelEditPart.getFigure().getBounds());
            targetBounds.x = (int) (klabel.getX() * scale);
            targetBounds.y = (int) (klabel.getY() * scale);
            
            ConnectionEditPart connectionEditPart = (ConnectionEditPart) labelEditPart.getParent();
            PointList bendPoints = getBendPoints((ElkEdge) parent, connectionEditPart.getFigure(), scale);
            EObject modelElement = connectionEditPart.getNotationView().getElement();
            EdgeLabelPlacement labelPlacement = klabel.getProperty(CoreOptions.EDGE_LABELS_PLACEMENT);
            
            // for labels of the opposite reference of an ecore reference,
            // the list of bend points must be reversed
            if (modelElement instanceof EReference && labelPlacement == EdgeLabelPlacement.TAIL) {
                bendPoints = bendPoints.getCopy();
                bendPoints.reverse();
            }
            
            // get the referencePoint for the label
            int fromEnd, keyPoint = ConnectionLocator.MIDDLE;
            if (labelEditPart instanceof LabelEditPart) {
                keyPoint = ((LabelEditPart) labelEditPart).getKeyPoint();
            }
            switch (keyPoint) {
            case ConnectionLocator.SOURCE:
                fromEnd = SOURCE_LOCATION;
                break;
            case ConnectionLocator.TARGET:
                fromEnd = TARGET_LOCATION;
                break;
            default:
                fromEnd = MIDDLE_LOCATION;
                break;
            }
            Point refPoint = PointListUtilities.calculatePointRelativeToLine(bendPoints, 0, fromEnd, true);
            
            // get the new relative location
            Point normalPoint = offsetFromRelativeCoordinate(targetBounds, bendPoints, refPoint);
            if (normalPoint != null) {
                command.addShapeLayout((View) labelEditPart.getModel(), normalPoint, null);
            }
        }
    }

    /**
     * Transform the bend points of the given edge layout into a point list, reusing existing ones
     * if possible. The source and target points of the edge layout are included in the point list.
     * 
     * @param edge
     *            the edge for which to fetch bend points
     * @param isSplineEdge
     *            indicates whether the connection supports splines
     * @return point list with the bend points of the edge layout
     * @param scale
     *            scale factor for coordinates
     */
    private PointList getBendPoints(final ElkEdge edge, final IFigure edgeFigure, final double scale) {
        // This assumes that the edge has at least one edge section, which by this point it should
        ElkEdgeSection edgeSection = edge.getSections().get(0);
        PointList pointList = pointListMap.get(edgeSection);
        if (pointList == null) {
            KVectorChain bendPoints = ElkUtil.createVectorChain(edgeSection);

            // for connections that support splines the control points are passed without change
            boolean approx = handleSplineConnection(edgeFigure, edge.getProperty(CoreOptions.EDGE_ROUTING));
            
            // in other cases an approximation is used
            if (approx && bendPoints.size() >= 1) {
                bendPoints = ElkMath.approximateBezierSpline(bendPoints);
            }

            bendPoints.scale(scale);
            pointList = new PointList(bendPoints.size() + 2);
            for (KVector bendPoint : bendPoints) {
                pointList.addPoint((int) bendPoint.x, (int) bendPoint.y);
            }

            pointListMap.put(edgeSection, pointList);
        }
        return pointList;
    }

    /**
     * class name of the ELK SplineConnection.
     * 
     * TODO: This class doesn't exist anymore...
     */
    private static final String SPLINE_CONNECTION
            = "org.eclipse.elk.core.model.gmf.figures.SplineConnection";

    /**
     * Handle the ELK SplineConnection class without a direct reference to it. Reflection is used
     * to avoid a dependency to its containing plugin.
     * 
     * @param edgeFigure
     *            the edge figure instance
     * @param edgeRouting
     *            the edge routing returned by the layout algorithm
     * @return {@code true} if an approximation should be used to represent the spline
     */
    private static boolean handleSplineConnection(final IFigure edgeFigure,
            final EdgeRouting edgeRouting) {
        boolean isSC;
        Class<?> clazz = edgeFigure.getClass();
        do {
            String canonicalName = clazz.getCanonicalName();
            // in some cases, eg anonymous classes, the canonical name may be null
            isSC = canonicalName != null && canonicalName.equals(SPLINE_CONNECTION);
            clazz = clazz.getSuperclass();
        } while (!isSC && clazz != null);
        if (isSC) {
            clazz = edgeFigure.getClass();
            try {
                if (edgeRouting == EdgeRouting.SPLINES) {
                    // SplineConnection.SPLINE_CUBIC
                    clazz.getMethod("setSplineMode", int.class).invoke(edgeFigure, 1);
                } else {
                    // SplineConnection.SPLINE_OFF
                    clazz.getMethod("setSplineMode", int.class).invoke(edgeFigure, 0);
                }
                return false;
            } catch (Exception exception) {
                throw new WrappedException(exception);
            }
        }
        // no spline connection class, but spline representation is requested
        return edgeRouting == EdgeRouting.SPLINES;
    }

    /**
     * <!-- CHECKSTYLEOFF LineLength --> Calculates the label offset from the reference point given
     * the label bounds and a points list. This code has been copied and adapted from
     * {@link org.eclipse.gmf.runtime.diagram.ui.internal.figures.LabelHelper#offsetFromRelativeCoordinate(IFigure, Rectangle, PointList, Point)}
     * ,
     * {@link org.eclipse.gmf.runtime.diagram.ui.internal.figures.LabelHelper#normalizeRelativePointToPointOnLine(PointList, Point, Point)}
     * , and
     * {@link org.eclipse.gmf.runtime.diagram.ui.internal.figures.LabelHelper#getOrthogonalDistances(LineSeg, Point, Point)}
     * .
     * 
     * <!-- CHECKSTYLEON LineLength -->
     * 
     * @param bounds
     *            the {@code Rectangle} that is the bounding box of the label
     * @param points
     *            the {@code PointList} that the label offset is relative to
     * @param therefPoint
     *            the {@code Point} that is the reference point that the offset is based on, or
     *            {@code null}
     * @return a {@code Point} which represents a value offset from the {@code refPoint} point
     *         oriented based on the nearest line segment, or {@code null} if no such point can be
     *         determined
     */
    @SuppressWarnings("restriction")
    public static Point offsetFromRelativeCoordinate(final Rectangle bounds,
            final PointList points, final Point therefPoint) {
        Point refPoint = therefPoint;
        if (refPoint == null) {
            refPoint = points.getFirstPoint();
        }
        // compensate for the fact that we are using the figure center
        bounds.translate(bounds.width / 2, bounds.height / 2);
        Point offset = new Point(bounds.x - refPoint.x, bounds.y - refPoint.y);
        // calculate slope of line
        if (points.size() == 1) {
            // this is a node...
            return offset;
        } else if (points.size() >= 2) {
            // this is an edge...
            int segIndex = PointListUtilities.findNearestLineSegIndexOfPoint(points, refPoint);
            @SuppressWarnings("rawtypes")
            List segmentsList = PointListUtilities.getLineSegments(points);
            if (segIndex <= 0) {
                segIndex = 0;
            } else if (segIndex > segmentsList.size()) {
                segIndex = segmentsList.size() - 1;
            } else {
                segIndex--;
            }
            LineSeg segment = (LineSeg) segmentsList.get(segIndex);
            Point normalOffset = null;
            if (segment.isHorizontal()) {
                if (segment.getOrigin().x > segment.getTerminus().x) {
                    normalOffset = offset.getNegated();
                    return normalOffset;
                } else {
                    normalOffset = offset;
                    return normalOffset;
                }
            } else if (segment.isVertical()) {
                if (segment.getOrigin().y < segment.getTerminus().y) {
                    normalOffset = offset.scale(-1, 1).transpose();
                    return normalOffset;
                } else {
                    normalOffset = offset.scale(1, -1).transpose();
                    return normalOffset;
                }
            } else {
                Point offsetRefPoint = refPoint.getTranslated(offset);
                LineSeg parallelSeg = segment.getParallelLineSegThroughPoint(offsetRefPoint);
                Point p1 = parallelSeg.perpIntersect(refPoint.x, refPoint.y);
                double dx = p1.getDistance(offsetRefPoint) * ((p1.x > offsetRefPoint.x) ? -1 : 1);
                double dy = p1.getDistance(refPoint) * ((p1.y < refPoint.y) ? -1 : 1);
                Point orth = new PrecisionPoint(dx, dy);
                // reflection in the y axis
                if (segment.getOrigin().x > segment.getTerminus().x) {
                    orth = orth.scale(-1, -1);
                }
                return orth;
            }
        }
        return null;
    }

}
