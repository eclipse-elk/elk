/**
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 */
package org.eclipse.elk.graph.util;

import java.util.Map;

import org.eclipse.elk.graph.*;

import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;

import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.ecore.EPackage;

import org.eclipse.emf.ecore.util.Switch;

/**
 * <!-- begin-user-doc -->
 * The <b>Switch</b> for the model's inheritance hierarchy.
 * It supports the call {@link #doSwitch(EObject) doSwitch(object)}
 * to invoke the <code>caseXXX</code> method for each class of the model,
 * starting with the actual class of the object
 * and proceeding up the inheritance hierarchy
 * until a non-null result is returned,
 * which is the result of the switch.
 * <!-- end-user-doc -->
 * @see org.eclipse.elk.graph.ElkGraphPackage
 * @generated
 */
public class ElkGraphSwitch<T> extends Switch<T> {
    /**
     * The cached model package
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    protected static ElkGraphPackage modelPackage;

    /**
     * Creates an instance of the switch.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkGraphSwitch() {
        if (modelPackage == null) {
            modelPackage = ElkGraphPackage.eINSTANCE;
        }
    }

    /**
     * Checks whether this is a switch for the given package.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @param ePackage the package in question.
     * @return whether this is a switch for the given package.
     * @generated
     */
    @Override
    protected boolean isSwitchFor(EPackage ePackage) {
        return ePackage == modelPackage;
    }

    /**
     * Calls <code>caseXXX</code> for each class of the model until one returns a non null result; it yields that result.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @return the first non-null result returned by a <code>caseXXX</code> call.
     * @generated
     */
    @Override
    protected T doSwitch(int classifierID, EObject theEObject) {
        switch (classifierID) {
            case ElkGraphPackage.IPROPERTY_HOLDER: {
                IPropertyHolder iPropertyHolder = (IPropertyHolder)theEObject;
                T result = caseIPropertyHolder(iPropertyHolder);
                if (result == null) result = defaultCase(theEObject);
                return result;
            }
            case ElkGraphPackage.EMAP_PROPERTY_HOLDER: {
                EMapPropertyHolder eMapPropertyHolder = (EMapPropertyHolder)theEObject;
                T result = caseEMapPropertyHolder(eMapPropertyHolder);
                if (result == null) result = caseIPropertyHolder(eMapPropertyHolder);
                if (result == null) result = defaultCase(theEObject);
                return result;
            }
            case ElkGraphPackage.ELK_GRAPH_ELEMENT: {
                ElkGraphElement elkGraphElement = (ElkGraphElement)theEObject;
                T result = caseElkGraphElement(elkGraphElement);
                if (result == null) result = caseEMapPropertyHolder(elkGraphElement);
                if (result == null) result = caseIPropertyHolder(elkGraphElement);
                if (result == null) result = defaultCase(theEObject);
                return result;
            }
            case ElkGraphPackage.ELK_SHAPE: {
                ElkShape elkShape = (ElkShape)theEObject;
                T result = caseElkShape(elkShape);
                if (result == null) result = caseElkGraphElement(elkShape);
                if (result == null) result = caseEMapPropertyHolder(elkShape);
                if (result == null) result = caseIPropertyHolder(elkShape);
                if (result == null) result = defaultCase(theEObject);
                return result;
            }
            case ElkGraphPackage.ELK_LABEL: {
                ElkLabel elkLabel = (ElkLabel)theEObject;
                T result = caseElkLabel(elkLabel);
                if (result == null) result = caseElkShape(elkLabel);
                if (result == null) result = caseElkGraphElement(elkLabel);
                if (result == null) result = caseEMapPropertyHolder(elkLabel);
                if (result == null) result = caseIPropertyHolder(elkLabel);
                if (result == null) result = defaultCase(theEObject);
                return result;
            }
            case ElkGraphPackage.ELK_CONNECTABLE_SHAPE: {
                ElkConnectableShape elkConnectableShape = (ElkConnectableShape)theEObject;
                T result = caseElkConnectableShape(elkConnectableShape);
                if (result == null) result = caseElkShape(elkConnectableShape);
                if (result == null) result = caseElkGraphElement(elkConnectableShape);
                if (result == null) result = caseEMapPropertyHolder(elkConnectableShape);
                if (result == null) result = caseIPropertyHolder(elkConnectableShape);
                if (result == null) result = defaultCase(theEObject);
                return result;
            }
            case ElkGraphPackage.ELK_NODE: {
                ElkNode elkNode = (ElkNode)theEObject;
                T result = caseElkNode(elkNode);
                if (result == null) result = caseElkConnectableShape(elkNode);
                if (result == null) result = caseElkShape(elkNode);
                if (result == null) result = caseElkGraphElement(elkNode);
                if (result == null) result = caseEMapPropertyHolder(elkNode);
                if (result == null) result = caseIPropertyHolder(elkNode);
                if (result == null) result = defaultCase(theEObject);
                return result;
            }
            case ElkGraphPackage.ELK_PORT: {
                ElkPort elkPort = (ElkPort)theEObject;
                T result = caseElkPort(elkPort);
                if (result == null) result = caseElkConnectableShape(elkPort);
                if (result == null) result = caseElkShape(elkPort);
                if (result == null) result = caseElkGraphElement(elkPort);
                if (result == null) result = caseEMapPropertyHolder(elkPort);
                if (result == null) result = caseIPropertyHolder(elkPort);
                if (result == null) result = defaultCase(theEObject);
                return result;
            }
            case ElkGraphPackage.ELK_EDGE: {
                ElkEdge elkEdge = (ElkEdge)theEObject;
                T result = caseElkEdge(elkEdge);
                if (result == null) result = caseElkGraphElement(elkEdge);
                if (result == null) result = caseEMapPropertyHolder(elkEdge);
                if (result == null) result = caseIPropertyHolder(elkEdge);
                if (result == null) result = defaultCase(theEObject);
                return result;
            }
            case ElkGraphPackage.ELK_BEND_POINT: {
                ElkBendPoint elkBendPoint = (ElkBendPoint)theEObject;
                T result = caseElkBendPoint(elkBendPoint);
                if (result == null) result = defaultCase(theEObject);
                return result;
            }
            case ElkGraphPackage.ELK_EDGE_SECTION: {
                ElkEdgeSection elkEdgeSection = (ElkEdgeSection)theEObject;
                T result = caseElkEdgeSection(elkEdgeSection);
                if (result == null) result = caseEMapPropertyHolder(elkEdgeSection);
                if (result == null) result = caseIPropertyHolder(elkEdgeSection);
                if (result == null) result = defaultCase(theEObject);
                return result;
            }
            case ElkGraphPackage.ELK_PROPERTY_TO_VALUE_MAP_ENTRY: {
                @SuppressWarnings("unchecked") Map.Entry<IProperty<?>, Object> elkPropertyToValueMapEntry = (Map.Entry<IProperty<?>, Object>)theEObject;
                T result = caseElkPropertyToValueMapEntry(elkPropertyToValueMapEntry);
                if (result == null) result = defaultCase(theEObject);
                return result;
            }
            default: return defaultCase(theEObject);
        }
    }

    /**
     * Returns the result of interpreting the object as an instance of '<em>IProperty Holder</em>'.
     * <!-- begin-user-doc -->
     * This implementation returns null;
     * returning a non-null result will terminate the switch.
     * <!-- end-user-doc -->
     * @param object the target of the switch.
     * @return the result of interpreting the object as an instance of '<em>IProperty Holder</em>'.
     * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
     * @generated
     */
    public T caseIPropertyHolder(IPropertyHolder object) {
        return null;
    }

    /**
     * Returns the result of interpreting the object as an instance of '<em>EMap Property Holder</em>'.
     * <!-- begin-user-doc -->
     * This implementation returns null;
     * returning a non-null result will terminate the switch.
     * <!-- end-user-doc -->
     * @param object the target of the switch.
     * @return the result of interpreting the object as an instance of '<em>EMap Property Holder</em>'.
     * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
     * @generated
     */
    public T caseEMapPropertyHolder(EMapPropertyHolder object) {
        return null;
    }

    /**
     * Returns the result of interpreting the object as an instance of '<em>Element</em>'.
     * <!-- begin-user-doc -->
     * This implementation returns null;
     * returning a non-null result will terminate the switch.
     * <!-- end-user-doc -->
     * @param object the target of the switch.
     * @return the result of interpreting the object as an instance of '<em>Element</em>'.
     * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
     * @generated
     */
    public T caseElkGraphElement(ElkGraphElement object) {
        return null;
    }

    /**
     * Returns the result of interpreting the object as an instance of '<em>Elk Shape</em>'.
     * <!-- begin-user-doc -->
     * This implementation returns null;
     * returning a non-null result will terminate the switch.
     * <!-- end-user-doc -->
     * @param object the target of the switch.
     * @return the result of interpreting the object as an instance of '<em>Elk Shape</em>'.
     * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
     * @generated
     */
    public T caseElkShape(ElkShape object) {
        return null;
    }

    /**
     * Returns the result of interpreting the object as an instance of '<em>Elk Label</em>'.
     * <!-- begin-user-doc -->
     * This implementation returns null;
     * returning a non-null result will terminate the switch.
     * <!-- end-user-doc -->
     * @param object the target of the switch.
     * @return the result of interpreting the object as an instance of '<em>Elk Label</em>'.
     * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
     * @generated
     */
    public T caseElkLabel(ElkLabel object) {
        return null;
    }

    /**
     * Returns the result of interpreting the object as an instance of '<em>Elk Connectable Shape</em>'.
     * <!-- begin-user-doc -->
     * This implementation returns null;
     * returning a non-null result will terminate the switch.
     * <!-- end-user-doc -->
     * @param object the target of the switch.
     * @return the result of interpreting the object as an instance of '<em>Elk Connectable Shape</em>'.
     * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
     * @generated
     */
    public T caseElkConnectableShape(ElkConnectableShape object) {
        return null;
    }

    /**
     * Returns the result of interpreting the object as an instance of '<em>Elk Node</em>'.
     * <!-- begin-user-doc -->
     * This implementation returns null;
     * returning a non-null result will terminate the switch.
     * <!-- end-user-doc -->
     * @param object the target of the switch.
     * @return the result of interpreting the object as an instance of '<em>Elk Node</em>'.
     * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
     * @generated
     */
    public T caseElkNode(ElkNode object) {
        return null;
    }

    /**
     * Returns the result of interpreting the object as an instance of '<em>Elk Port</em>'.
     * <!-- begin-user-doc -->
     * This implementation returns null;
     * returning a non-null result will terminate the switch.
     * <!-- end-user-doc -->
     * @param object the target of the switch.
     * @return the result of interpreting the object as an instance of '<em>Elk Port</em>'.
     * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
     * @generated
     */
    public T caseElkPort(ElkPort object) {
        return null;
    }

    /**
     * Returns the result of interpreting the object as an instance of '<em>Elk Edge</em>'.
     * <!-- begin-user-doc -->
     * This implementation returns null;
     * returning a non-null result will terminate the switch.
     * <!-- end-user-doc -->
     * @param object the target of the switch.
     * @return the result of interpreting the object as an instance of '<em>Elk Edge</em>'.
     * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
     * @generated
     */
    public T caseElkEdge(ElkEdge object) {
        return null;
    }

    /**
     * Returns the result of interpreting the object as an instance of '<em>Elk Bend Point</em>'.
     * <!-- begin-user-doc -->
     * This implementation returns null;
     * returning a non-null result will terminate the switch.
     * <!-- end-user-doc -->
     * @param object the target of the switch.
     * @return the result of interpreting the object as an instance of '<em>Elk Bend Point</em>'.
     * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
     * @generated
     */
    public T caseElkBendPoint(ElkBendPoint object) {
        return null;
    }

    /**
     * Returns the result of interpreting the object as an instance of '<em>Elk Edge Section</em>'.
     * <!-- begin-user-doc -->
     * This implementation returns null;
     * returning a non-null result will terminate the switch.
     * <!-- end-user-doc -->
     * @param object the target of the switch.
     * @return the result of interpreting the object as an instance of '<em>Elk Edge Section</em>'.
     * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
     * @generated
     */
    public T caseElkEdgeSection(ElkEdgeSection object) {
        return null;
    }

    /**
     * Returns the result of interpreting the object as an instance of '<em>Elk Property To Value Map Entry</em>'.
     * <!-- begin-user-doc -->
     * This implementation returns null;
     * returning a non-null result will terminate the switch.
     * <!-- end-user-doc -->
     * @param object the target of the switch.
     * @return the result of interpreting the object as an instance of '<em>Elk Property To Value Map Entry</em>'.
     * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
     * @generated
     */
    public T caseElkPropertyToValueMapEntry(Map.Entry<IProperty<?>, Object> object) {
        return null;
    }

    /**
     * Returns the result of interpreting the object as an instance of '<em>EObject</em>'.
     * <!-- begin-user-doc -->
     * This implementation returns null;
     * returning a non-null result will terminate the switch, but this is the last case anyway.
     * <!-- end-user-doc -->
     * @param object the target of the switch.
     * @return the result of interpreting the object as an instance of '<em>EObject</em>'.
     * @see #doSwitch(org.eclipse.emf.ecore.EObject)
     * @generated
     */
    @Override
    public T defaultCase(EObject object) {
        return null;
    }

} //ElkGraphSwitch
