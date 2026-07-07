/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later  
 *******************************************************************************/
package org.eclipse.elk.alg.topdownpacking;

import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;
import org.eclipse.emf.common.notify.Adapter;
import org.eclipse.emf.common.notify.Notification;
import org.eclipse.emf.common.util.EList;
import org.eclipse.emf.common.util.EMap;
import org.eclipse.emf.common.util.TreeIterator;
import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.ecore.EOperation;
import org.eclipse.emf.ecore.EReference;
import org.eclipse.emf.ecore.EStructuralFeature;
import org.eclipse.emf.ecore.resource.Resource;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * An {@link ElkNode} that stores its children in a two-dimensional grid.
 *
 */
public class GridElkNode implements ElkNode, Grid<ElkNode> {
    
    // The grid is stored as a list of rows and each row contains a list of column entries for 
    // that row.
    private List<List<ElkNode>> grid;
    private int rows = 0;
    private int cols = 0;
    
    // wrapped ElkNode
    private ElkNode node;

    /**
     * Wraps an existing {@link ElkNode} in a grid wrapper.
     * @param elkNode the base elk node
     */
    public GridElkNode(ElkNode elkNode) {
        this.node = elkNode;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void put(int col, int row, ElkNode item) throws IndexOutOfBoundsException {
        if (col >= cols || row >= rows) {
            throw new IndexOutOfBoundsException("The grid has a size of (" 
        + cols + " x " + rows + "). The requested index of (" + col + " x " + row + ") is out of bounds.");
        }
        
        grid.get(row).add(col, item);
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ElkNode get(int col, int row) throws IndexOutOfBoundsException {
        if (col >= cols || row >= rows) {
            throw new IndexOutOfBoundsException("The grid has a size of (" 
        + cols + " x " + rows + "). The requested index of (" + col + " x " + row + ") is out of bounds.");
        }
        
        return grid.get(row).get(col);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public List<ElkNode> getRow(int row) throws IndexOutOfBoundsException {
        if (row >= rows) {
            throw new IndexOutOfBoundsException("The grid has a size of (" 
        + cols + " x " + rows + "). The requested row " + row + ") is out of bounds.");
        }
        return grid.get(row);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public List<ElkNode> getColumn(int col) throws IndexOutOfBoundsException {
        if (col >= cols) {
            throw new IndexOutOfBoundsException("The grid has a size of (" 
        + cols + " x " + rows + "). The requested column " + col + ") is out of bounds.");
        }
        ArrayList<ElkNode> resultList = new ArrayList<>();
        for (int i = 0; i < rows; i++) {
            resultList.add(grid.get(i).get(col));
        }
        return resultList;
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getColumns() {
        return cols;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getRows() {
        return rows;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setGridSize(int cols, int rows) {
        // initialize the backing arraylists and fill them with null
        this.cols = cols;
        this.rows = rows;
        this.grid = new ArrayList<List<ElkNode>>(rows);
        
        for (int i = 0; i < rows; i++) {
            grid.add(new ArrayList<ElkNode>());
            for (int j = 0; j < cols; j++) {
                grid.get(i).add(null);
            }
        }
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EList<ElkEdge> getOutgoingEdges() {
        return node.getOutgoingEdges();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EList<ElkEdge> getIncomingEdges() {
        return node.getIncomingEdges();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getHeight() {
        return node.getHeight();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setHeight(double value) {
        node.setHeight(value);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getWidth() {
        return node.getWidth();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setWidth(double value) {
        node.setWidth(value);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getX() {
        return node.getX();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setX(double value) {
        node.setX(value);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getY() {
        return node.getY();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setY(double value) {
        node.setY(value);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setDimensions(double width, double height) {
        node.setDimensions(width, height);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setLocation(double x, double y) {
        node.setLocation(x, y);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EList<ElkLabel> getLabels() {
        return node.getLabels();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String getIdentifier() {
        return node.getIdentifier();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setIdentifier(String value) {
        node.setIdentifier(value);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EMap<IProperty<?>, Object> getProperties() {
        return node.getProperties();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EClass eClass() {
        return node.eClass();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Resource eResource() {
        return node.eResource();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EObject eContainer() {
        return node.eContainer();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EStructuralFeature eContainingFeature() {
        return node.eContainingFeature();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EReference eContainmentFeature() {
        return node.eContainmentFeature();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EList<EObject> eContents() {
        return node.eContents();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public TreeIterator<EObject> eAllContents() {
        return node.eAllContents();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean eIsProxy() {
        return node.eIsProxy();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EList<EObject> eCrossReferences() {
        return node.eCrossReferences();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Object eGet(EStructuralFeature feature) {
        return node.eGet(feature);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Object eGet(EStructuralFeature feature, boolean resolve) {
        return node.eGet(feature, resolve);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void eSet(EStructuralFeature feature, Object newValue) {
        node.eSet(feature, newValue);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean eIsSet(EStructuralFeature feature) {
        return node.eIsSet(feature);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void eUnset(EStructuralFeature feature) {
        node.eUnset(feature);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Object eInvoke(EOperation operation, EList<?> arguments)
            throws InvocationTargetException {
        return node.eInvoke(operation, arguments);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EList<Adapter> eAdapters() {
        return node.eAdapters();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean eDeliver() {
        return node.eDeliver();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void eSetDeliver(boolean deliver) {
        node.eSetDeliver(deliver);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void eNotify(Notification notification) {
        node.eNotify(notification);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public <T> IPropertyHolder setProperty(IProperty<? super T> property, T value) {
        return node.setProperty(property, value);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public <T> T getProperty(IProperty<T> property) {
        return node.getProperty(property);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasProperty(IProperty<?> property) {
        return node.hasProperty(property);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public IPropertyHolder copyProperties(IPropertyHolder holder) {
        return node.copyProperties(holder);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Map<IProperty<?>, Object> getAllProperties() {
        return node.getAllProperties();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EList<ElkPort> getPorts() {
        return node.getPorts();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EList<ElkNode> getChildren() {
        return node.getChildren();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ElkNode getParent() {
        return node.getParent();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setParent(ElkNode value) {
        node.setParent(value);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public EList<ElkEdge> getContainedEdges() {
        return node.getContainedEdges();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean isHierarchical() {
        return node.isHierarchical();
    }

}
